site_name: SLOTH 
repo_url: https://www-git-cad.intra.cea.fr/DEC/collaboratif/ci230846/DEV_PROJECT/sloth.git
use_directory_urls: false
markdown_extensions:  
  - toc:
      permalink: true
  - attr_list
  - md_in_html
  - pymdownx.highlight:
      anchor_linenums: true
      line_spans: __span
      pygments_lang_class: true
  - pymdownx.inlinehilite
  - pymdownx.snippets:
      base_path: $relative
  - pymdownx.superfences
  - def_list
  - pymdownx.tasklist:
      custom_checkbox: true
  - pymdownx.arithmatex:
      generic: true

extra_javascript:
  - javascripts/mathjax.js
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js

plugins:
  - offline
  - search
  - with-pdf

theme: 
  name: material
  logo: img/sloth_bl.png
  icon:
    repo: fontawesome/brands/gitlab
    previous: fontawesome/solid/angle-left
    next: fontawesome/solid/angle-right

  features:
    - navigation.tabs
    # - navigation.sections
    # - navigation.indexes
    - navigation.top
    - navigation.instant
    - navigation.instant.progress
    - navigation.tracking
    - content.code.copy
  palette:
    primary: black
    accent: 'red'
########################################
########################################
docs_dir: '@DOCS_DIR@'
########################################
########################################
nav:
  - Home: index.md
  - Getting Started:
    - Installing SLOTH: Started/installation.md
    - Creating a Phase-Field application: Started/howto.md
    - Examples: 
      - Synthesis: Started/examples.md
      @EXAMPLE_LIST@
  - Applications:
    - Oxygen Thermal Diffusion: Applications/otd.md
    - Incipient Melting: Applications/melting.md

  - Documentation: 
      - User Manual: 
        - Introduction: Documentation/User/user_manual.md
        - Model 2: Documentation/User/model2.md
        - Model 3: Documentation/User/model3.md
      - Modelling Description: 
        - Introduction: Documentation/Physical/physical_description.md
        - Model 1: Documentation/Physical/model1.md
        - Model 2: Documentation/Physical/model2.md
        - Model 3: Documentation/Physical/model3.md
      - Code Documentation: Documentation/Code/code_description.md
  - V&V: VerifValid/VV.md
  - References: References/Ref.md
  - About: about.md
# rajouter la recherche des models, des examples, des post-pro, pandox plugin