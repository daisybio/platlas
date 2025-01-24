#test for system dependencies
renv_lock_contents <- readLines("renv.lock")

# Extract the package names from the lock file
package_lines <- grep("Package", renv_lock_contents, value = TRUE)

# Extract the package names from the lines
package_names <- str_extract(package_lines, "(?<=\"Package\": \")[^\"]+")

# Print the package names
package_names <- package_names[!is.na(package_names)]

system_dependencies <- list()

# Iterate over the package names
#for (pkg in package_names) {
  # Get the package dependencies for the current package
pdb <- available.packages()
dependencies <- tools::package_dependencies(package_names,pdb,which = "Depends")
#library(sysreqs)
#sysdeps <- sysreqs::sysreqs("htmltools")
  # Extract the system dependencies
  #sysdeps <- dependencies$System

for (pkg in package_names) {
  #print(pkg)
  query <- maketools::package_sysdeps(pkg)
  sysdep <- query$shlib
  system_dependencies[[pkg]] <- sysdep
}

sys_deps <- unlist(system_dependencies)
sys_deps <- unique(sys_deps)
  # Add the system dependencies to the list
  #system_dependencies[[pkg]] <- sysdeps
#}
#linux_sysdeps <-  maketools::sysdep_map()
# if ("zlib" %in% SystemRequirements()) {
#   zlib_libname <- SystemRequirements()[["zlib"]]
#   print(paste("The zlib library is available and the corresponding library name is", zlib_libname))
# } else {
#   print("The zlib library is not available on your system")
# }

# Print the system dependencies
#system_dependencies