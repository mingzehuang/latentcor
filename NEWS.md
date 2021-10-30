# latentcor 1.1.0

* Rename main function from `estR` to `latentcor`.
* Add help function `get_types` to automatically generate types for data matrix.
* Add more documentation on parameters.
* Add more documentation for speed comparison between original method and approximation.

# latentcor 1.2.0

* Add user-defined `use.nearPD`, so that user can decide if latent correlation matrix should be adjusted to be positive definite automatically.
* Remove redundant code for checking positive definiteness.
* Minor correction on types detection to accommodate `NA` values. 
