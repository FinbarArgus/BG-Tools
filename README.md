# Generate CellML from BG-RDF models

> [!WARNING]
> These tools are actively being developed and nothing should be considered stable until version `1.0` has been released.

## Installation

### As a package

```sh
uv add https://github.com/CellDL/BG-Tools/releases/download/v0.9.2/bg2cellml-0.9.2-py3-none-any.whl
```

### Development

```sh
$ git clone https://github.com/CellDL/BG-Tools
$ cd BG-Tools
$ git submodule update --init --recursive
$ uv sync --all-groups
```

## Usage

See `cli/bgtool/__main__.py` for details.

### A module for libcuflynx and PhLynx

```sh
$ PYTHONPATH=cli python -m bgtool --module glycogen/bondgraph_v01 --output OUT_DIR glycogen.celldl.svg
```

writes a module version as a module library (such as circulatory-autogen-modules) holds one, in the formats
of libcuflynx's `schemas/module_config.schema.json`, PhLynx's keys: `<type>_<subtype>_modules.cellml`,
`_units.cellml`, `_modules_config.json` and `parameterisations/default/default_parameters.csv`, every value
of the bond graph a parameter there. Each chemostat species and compartment (all its copies) is a general
port, `chemostat_<species>_<location>`, of its potential (a boundary condition) and the module's uptake of
the species. See `bg2cellml/module_export.py`. `tests/test_module_export.py` builds the module with
libcuflynx and runs it against bg2cellml's own CellML when `LIBCUFLYNX_SRC` (circulatory_autogen's `src/`)
and `LIBCUFLYNX_PYTHON` (a Python with its dependencies) are set.

### A variable mapping for libcuflynx's methods LaTeX

```sh
$ PYTHONPATH=cli python -m bgtool --variable-mapping --output OUT_DIR glycogen.celldl.svg
```

also writes `OUT_DIR/<name>_variable_mapping.csv` beside the CellML, and `--module` writes
`<type>_<subtype>_variable_mapping.csv` with the module. It is libcuflynx's variable mapping
(`cuflynx-methods-latex --mapping`): `variable_name` (`component/variable`), `latex`, `kind`, `component` and
`units`, with one more column, `bondgraph`, each variable's bond graph name (`q_G6P_c`). A variable keeping its
bond graph name is written in bond graph notation (`q^{\mathrm{G6P}}_{c}`); one the model names its own way
(bgf:NamingRule, bgf:VariableName) by libcuflynx's default rule for its new name. See
`bg2cellml/variable_mapping.py`.

## Example models:

* [An RC electrical circuit](https://github.com/CellDL/BG-RDF/blob/main/examples/example_RC.ttl) 
* [A simple biochemical reaction](https://github.com/CellDL/BG-RDF/blob/main/examples/example_B1.ttl)
* [All examples](https://github.com/CellDL/BG-RDF/blob/main/examples)

### Examples in a Jupyter notebook:

1. Clone the repository and its sub-modules as above.
2. Install with `uv sync --group notebook`
3. Start JupyterLab with `uv run jupyter lab`
4. Open `./notebooks/bg2cellml.ipynb` in Jupyter.
