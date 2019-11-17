# EPICO-data-loading-scripts install instructions

In order to use the available scripts you need to locally install the dependencies.
Either [cpanm](https://metacpan.org/pod/App::cpanminus) or [cpm](https://metacpan.org/pod/App::cpm) can used for that task.

## Local install
As official releases of these tools are available in [BSC INB DarkPAN](https://gitlab.bsc.es/inb/dartpan), you only need to run next command:

```bash
cpanm --mirror-only --mirror https://gitlab.bsc.es/inb/darkpan/raw/master/ --mirror https://cpan.metacpan.org/ BP::DCCLoader
```

If you prefer to clone this repo, you only have to run next command in order to have all the needed dependencies installed at `local`:

```bash
cpm install --resolver 02packages,https://gitlab.bsc.es/inb/darkpan/raw/master/ --resolver metadb
```

## Install in a bigger project

If you project is managing its dependences using a `cpanfile` and `cpm`, the `cpanfile` should contain next line:

```cpanfile
requires 'BP::DCCLoader', 'v1.0.2';
```

and use either `cpm` or `cpanm` to install all the dependencies, using the
`--resolver` or `--mirror` flags to also search on BSC INB DarkPAN (as it was showed above).

