Here is the source code of paper "Negative Keyword-Aware Size-Constrained Community Search"

File folder "ciao_epinions_process" processes the source data of ciao and epinions.

File folder "citations_process" processes the source data of citaions.

File folder "dblp_process" processes the source data of dblp.

Files large_graph_gen.py and large_attribute_gen.py are the code for generating X * |V| large-scale graphs

Example of compiling and running (Linux)
Before running, you should:

Prepare the dataset and query files according to the templates
Check file names and file paths.\\
Then run the following commands,
```shell
$ make
$ main ciao v 0 30
# For running queries, five parameters should be specified:
# "bk" is the name of dataset
# "v" means that the experinment for negative size, "p" means that the experinment for positive query keyword, "n" means that the experinment for negative query keyword, "k" means that the experinment for k
# "0" means that we are running queries on varying algorithm
# "30" is the limition of runtime
```
