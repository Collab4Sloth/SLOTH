# SLOTH unit tests

Unit tests of the SLOTH kernel, written with [Catch2 v3](https://github.com/catchorg/Catch2).
They check one class or one function at a time: values against formulas
computed by hand, Jacobians against finite differences, error paths, and a few
small complete problems with exact solutions. The other directories of `tests/`
contain the calculation tests (complete simulations compared with reference
results).

## Dependencies

| Dependency | Needed for | Without it |
|---|---|---|
| Catch2 v3 | building the unit tests | the unit tests are skipped (CMake message), SLOTH and the other tests are built as usual |
| MFEM built with `MFEM_USE_EXCEPTIONS` | checking error paths (`SLOTH_CHECK_MFEM_ERROR`) | these checks are skipped with a warning, the other checks run |

The spack environment of SLOTH (`spack.yaml`) provides both
(`mfem +exceptions`, `catch2@3`). On macOS: `brew install catch2` (the MFEM of
Homebrew has no exceptions). Another Catch2 installation:
`-DCatch2_DIR=<prefix>/lib/cmake/Catch2`.

## Build and run

```bash
make unit_tests                         # unit tests only (make tests builds them too)
ctest -L Unit                           # all unit tests
ctest -L Integrators                    # one module
ctest -L Unit --output-on-failure

# An executable can be run directly, with the Catch2 options
mpirun -np 1 tests/Unit/Unit_Parameters                         # everything
mpirun -np 1 tests/Unit/Unit_Integrators "[residual]"           # tests with a tag
mpirun -np 1 tests/Unit/Unit_Integrators "AllenCahn*"           # tests by name
mpirun -np 1 tests/Unit/Unit_Parameters --list-tests
mpirun -np 1 tests/Unit/Unit_Parameters --success               # report successes too
```

Each module is run on 1 and 2 MPI ranks by default (`Unit.<Module>.np1`,
`Unit.<Module>.np2`); the CI runs at most 2 ranks.

## Layout

```
tests/Unit/
├── CMakeLists.txt           sloth_add_unit_test() and the list of modules
├── framework/
│   ├── sloth_unit_main.cpp  common main (MPI, Hypre, MFEM errors as exceptions, kernel timers)
│   ├── sloth_unit.hpp       helpers (see below)
│   └── test_mpi.cpp         checks of the harness itself
└── <Module>/                one directory per kernel module (names of kernel/include)
    └── test_<Class>.cpp     one file per class, or per group of classes
```

One executable per module (`Unit_<Module>`), built from all the files of the
module.

## Adding tests

**A file in an existing module**: create `<Module>/test_<Class>.cpp` and add it
to the `SOURCES` of the module in `CMakeLists.txt`.

**A new module**:

```cmake
sloth_add_unit_test(NAME Time
                    SOURCES Time/test_TimeDiscretization.cpp
                    NPROCS 1 2                                # optional, default 1 2
                    COEFFICIENTS_JSON Time/Coefficients.json) # optional, see below
```

`COEFFICIENTS_JSON`: coefficients generated at configure time by
`scripts/GenerateCoefficient.py` from a JSON file, as in a SLOTH calculation;
the generated headers are included by name (`#include "MobilityPhi.hpp"`, see
`Integrators/MobilityCoefficients.json`). Rerun `cmake .` after a change of the
JSON file.

**A test file**:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Parameters/Parameters.hpp"  // the unit under test
#include "sloth_unit.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("Parameters: get_param_value_or_default", "[Parameters]") {
  const Parameters params(Parameter("nb_iter", 10), Parameter("epsilon", 1.e-3));
  CHECK(params.get_param_value_or_default<int>("missing", 3) == 3);
  CHECK_THAT(params.get_param_value<double>("epsilon"), WithinRel(1.e-3, 1.e-12));

  SECTION("wrong type is an error") {
    SLOTH_CHECK_MFEM_ERROR_WITH(params.get_param_value<double>("nb_iter"),
                                "Invalid conversion");
  }
}
```

Tests that need complete SLOTH objects (`SPA`, `BCS`, `VAR`, operators,
problems...) include `"Sloth/sloth.hpp"` and `"Sloth/tests.hpp"` and use the
aliases of `Sloth2D` (`using namespace Sloth2D;`).

Useful Catch2 features: `CHECK` (goes on after a failure), `REQUIRE` (stops the
test case), `SECTION` (each section reruns the test case from the start, so
that the setup is written once), `TEST_CASE_METHOD` (fixture), `GENERATE` (one
run per value), `TEMPLATE_TEST_CASE` (one run per type), `CAPTURE`/`INFO`
(values printed with a failure), tag `[!mayfail]` (failure reported without
failing the executable: known defects).

## Helpers (`framework/sloth_unit.hpp`)

| Helper | Use |
|---|---|
| `SLOTH_CHECK_MFEM_ERROR(expr)` | `expr` ends in `mfem_error` (SLOTH errors, `MFEM_VERIFY`, `MFEM_ABORT`) |
| `SLOTH_CHECK_MFEM_ERROR_WITH(expr, "text")` | same, and the message contains `text` |
| `sloth_unit::point({x, y})` | `mfem::Vector` from a list of values |
| `sloth_unit::CoutCapture`, `CerrCapture` | what is written on `std::cout` / `std::cerr` while the object is alive |
| `sloth_unit::TemporaryFile(name, content)` | file in the temporary directory, removed at the end of the test |
| `sloth_unit::TemporaryDirectory(name)` | directory shared by all the ranks (collective), removed at the end of the test |
| `sloth_unit::UnitSquare(n, order)` | unit square, n x n quadrangles, H1 space of the given order, on all the ranks; usable as a fixture |

In the macros, an argument containing commas outside parentheses (braced
lists, template arguments) must be put in parentheses or in a variable first.

## Writing rules

- Expected values come from formulas or hand calculations, never from a
  previous SLOTH run.
- Fields exactly represented by the finite elements (affine fields in Q1,
  uniform fields) make the expected values exact, with no discretisation error.
- Jacobians are checked against centred finite differences of the residual.
- A test that needs a file writes it in a `TemporaryFile`/`TemporaryDirectory`:
  no file is left in the source or build tree.
- Comments describe what is tested, in the present tense; the history of a
  change goes in the commit message.
- MPI: every test runs on every rank, with the same order on all the ranks (the
  random seed of Catch2 is broadcast). Collective operations (`ParMesh`,
  `ComputeL2Error`, `MPI_Allreduce`...) are therefore safe; a check made on rank
  0 only (files written by rank 0) is written `if (mfem::Mpi::Root())`.

## Coverage

In a `Coverage` build (`envSloth.sh --coverage`):

```bash
make cc_reset                       # reset the coverage counters
ctest -L Unit                       # tests to measure (all of them: ctest)
make cc                             # report
python3 ../scripts/uncovered_functions.py Coverage/coverage.json -m Integrators
```

`make cc` writes `Coverage/coverage.html` (detailed, file by file),
`coverage.xml` and `coverage.json`, and prints the number of uncovered kernel
functions per module. The instantiations of a template are counted as one
function. Only the tests that have been run are taken into account. On GitHub
Actions (Linux job), the same table is written in the job summary.
