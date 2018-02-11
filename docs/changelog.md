# Changes in PMC tracking code

## v0.1 (2018-02-11)
*not a release!*

### Major changes
* Code structure is rewritten so that only one time slice is processed at each iteration. This allows for longer runs, because there is no need to store all time slices in memory.
* Due to this new logic, it is practically unfeasible to check if a vortex track has duration greater than a given threshold. This is now should be done at the postprocessing stage.
* File input is now done via netCDF interface
* Use `datetime` module
* Settings are stored in a `settings.conf` file instead of the main program
* Constants are stored in modules - thus reducing number of subroutine arguments

### Minor changes
* Add option to merge vortex depending on their respective sizes (#173f217)
* Add `great_circle()` function and trigonometric functions with arguments in degrees
* Remove (most) of the unused variables
* Fix certain array index bugs
