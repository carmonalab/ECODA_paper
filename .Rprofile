if (nzchar(Sys.getenv("PIXI_PROJECT_ROOT"))) {
  pixi_python <- Sys.which("python")
  
  if (nzchar(pixi_python)) {
    Sys.setenv(RETICULATE_PYTHON = pixi_python)
  }
} else {
  # Fallbacks for opening IDEs (RStudio/Positron) outside the active pixi shell
  if (.Platform$OS.type == "windows") {
    # Point directly to the default environment on Windows
    Sys.setenv(RETICULATE_PYTHON = file.path(getwd(), ".pixi/envs/default/Scripts/python.exe"))
  } else {
    # Point directly to the default environment on Mac/Linux HPC
    Sys.setenv(RETICULATE_PYTHON = file.path(getwd(), ".pixi/envs/default/bin/python"))
  }
}