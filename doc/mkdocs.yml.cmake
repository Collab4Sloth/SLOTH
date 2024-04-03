site_name: SLOTH 
repo_url: https://github.com/CINTROINI/sloth.git
use_directory_urls: false
markdown_extensions:  
  - toc:
      permalink: true
  - attr_list
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

plugins:
  - search
  - with-pdf

theme: 
  name: material
  logo: img/sloth.ico
  icon:
    repo: fontawesome/brands/github-alt
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
    primary: indigo
    accent: blue
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