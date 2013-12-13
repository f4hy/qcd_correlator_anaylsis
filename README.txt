Analysis code for analyzing Operators and Correlators from Lattice QCD
measurements

Explaination of each executable:


* alt_zfactor.py

  * Computes Z-factors (operator-overlaps) using fits to the
    diagonalized correlators
* binarywriter.py
  * Writes correlator data in a binary format to be used with older
    c++ analysis code
* build_corr.py
  * Convenience methods for creating correlator files given data files
* compare_cormatrix.py
  * Computes the chi-sqr differences of two entire correlator matricies
* configtimeobj.py
  * Class to store "config,time" indexed data such as correlator and
    operator data
* correlator.py
  * Extends the configtimeobj to have correlator specific methods,
    namely vevs
* determine_operators.py
  * Given a list of files for a correlator matrix, determine the list
    of operators from the filenames
* eigenvalues.py
  * Methods for reading Laph eigenvalues for glueball measurements
* fakecor_matrix.py
  * Methods to generate fake correlator maticies for testing
* fakecor.py
  * Method to generate fake correlator data
* fitfunctions.py
  * Functional forms to fit correlators to
* fit.py
  * Perform a covariant least squares fit to a correlator
* format_fit_results.py
  * Read in files for many fits and consolidate into a single file
* histo.py
  * Compute a histogram of a measurement for checking outliers/shape
* irreps.py
  * Figure out lattice irreps from continuum particle name
* jackknife.py
  * Basic jackknife functions for correlators (can probably be merged
    with correlator)
* level_identifier.py
  * Plot z-factor/rotation coeffs
* main.py
  * Compute effective masses. Should be renamed effective mass
* operator_tranlator.py
  * Translate operator names for printing in plots
* pandas_reader.py
  * read data files into pandas module objects
* particle_operators.py
  * Translating PDG names
* particles.py
  * List of particle names
* plot_files.py
  * Plot data files
* plot.py
  * Older plot library writing gnuplot files
* read_config_time_file.py
  * Older data file reader
* readinput.py
  * Simply library for reading various options from the user at
    command line
* test_cfgtimeobj.py
  * nit tests for cfgtimeobj
* tmin.py
  * Make tmin plots from a full correlator file
* vev.py
  * object for vevs to go with correlator objects
* zfactor.py
  * Compute zfactors(overlaps) using rotation coeffs, original
    correlator matrix
