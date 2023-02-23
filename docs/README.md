# PySEP Documentation

The PySEP documentation is built with Sphinx and ReadTheDocs. The official 
documentation can be found at https://adjtomo.github.io/pysep

## Building Docs Locally

In order to build the docs locally, you will first need to create a separate 
Conda environment with a few packages, you can do this by running:

``` bash
conda env create --file environment.yaml
conda activate pysep-docs
```

You can then run the make command to generate the .html files. You can find your 
local docs in the *_build/html* directory

```bash
make html
```

See your locally built documentation by opening *_build/html/index.html*

## Publishing Documentation

PySEP documentation is hosted on GitHub pages and currently must be built and
pushed manually (https://github.com/adjtomo/pysep/issues/74). Docs are built
from the ``gh-pages`` branch via GitHub actions.

If changes are made to the documentation, follow these instructions to push
changes to the live website.

1. Build the documentation locally following above instructions
2. Push documentation changes to ``gh-pages`` branch

	```bash
	mkdir /tmp/pysep-docs
	mv _build/html/* /tmp/pysep-docs
	cd ..
	git checkout gh-pages
	rm -rf *  # remove any existing files
	mv /tmp/pysep-docs/* .
	git add -A
	git commit -m # your commit message here
	git push --force origin gh-pages  # overwrite any existing fiels

3. Trigger the build action from https://github.com/adjtomo/pysep/actions. Look 
   for action 'pages build and deployment`. If not triggered automatically,
   click the button ``Re-run all jobs``
4. Check that your updates have been published

