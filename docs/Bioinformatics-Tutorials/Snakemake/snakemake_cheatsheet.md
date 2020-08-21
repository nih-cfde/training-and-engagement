# snakemake command cheatsheet

Snakemake | description
--- | ---
`snakemake -h` | help page with all Snakemake flags
`snakemake -n <rule name>` | a dry run of the rule commands without actually running them, good for testing
`snakemake -p <rule name>` | print shell commands to terminal as you run through the specified rule commands
`snakemake --delete-all-output <rule name>` | delete all the outputs up to a certain rule if specified, otherwise only output of first rule by default
