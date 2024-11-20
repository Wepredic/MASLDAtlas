library(reticulate)

reticulate::install_python(version = '3.9')


virtualenv_create(envname = "fibrosis_shiny",
                  version  = "3.9")
conda_list()
virtualenv_list()



use_virtualenv("fibrosis_shiny")

py_install("scanpy")
py_install("python-igraph")
py_install("leidenalg")
py_install("decoupler")
py_install("omnipath")
py_install("marsilea")
py_install("pydeseq2")
py_install("adjustText")
py_install("psutil")

