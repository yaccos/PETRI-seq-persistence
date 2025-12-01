Plan: Align Snakemake Layout

Adopt the Snakemake cookiecutter structure by moving configuration, environments, scripts, documentation, and outputs into the expected top-level folders. Update paths inside `workflow/Snakefile` and rule modules so rules reference the new locations, ensuring logs/results are separated from source and large resources sit under `resources/`.

Steps
1. Introduce `config/` and relocate `workflow/config.yaml`, `workflow/config_shiny.R`, then adjust `workflow/Snakefile` `configfile:` and any `include:` pointers.
2. Unlike the cookie-cutter Snakemake workflow, keep a global Conda environment for the entire workflow and store the `environment.yml` in the root of the repository.
3. Move user scripts into `workflow/scripts/` (Python/R) and clean stray `__pycache__` folders; fix rule `script:` paths accordingly.
4. Extract documentation and reports into `workflow/report/` (and `workflow/notebooks/` if notebooks exist or arrive later).
5. Migrate generated outputs and logs from `workflow/data/` and `workflow/results/` into top-level `results/` and `results/logs/`, updating rule `output` and `log` targets.
6. Consolidate reference genomes and raw FASTQs into `resources/` (or `resources/raw/`), updating any `input` directives that read from their old paths.

Further Considerations
1. Keep the folder `benchmark/` as-is. It is used for another branch of the repo and is under `.gitignore` and hence not part of the workflow.
