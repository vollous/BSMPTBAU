The required parameter lastline is not set.
CalcBAU calculates the baryon assymetry of the Universe 
it is called by

	./bin/CalcBAU model input output firstline lastline

or with arguments

	./bin/CalcBAU[arguments]

with the following arguments, ([*] are required arguments, others are optional):

argument                                default   description
--help                                            shows this menu
--model=                                          [*] model name
--input=                                          [*] input file (in tsv format)
--output=                                         [*] output file (in tsv format)
--firstline=                                      [*] line number of first line in input file
                                                      (expects line 1 to be a legend)
--lastline=                                       [*] line number of last line in input file
--thigh=                                300       high temperature [GeV]
--multistepmode=                        default   multi-step PT mode
                                                  default: default mode
                                                  0: single-step PT mode
                                                  >0 for multi-step PT modes:
                                                  1: tracing coverage
                                                  2: global minimum tracing coverage
                                                  auto: automatic mode
--num_pts=                              10        intermediate grid-size for default mode
--vwall=                                0.5       wall velocity: >0 user defined
                                                  -1: approximation
                                                  -2: upper bound
--perc_prbl=                            0.71      false vacuum fraction for percolation
--compl_prbl=                           0.01      false vacuum fraction for completion
--trans_temp=                           perc      transition temperature, options are:
                                                  crit: critical temperature
                                                  nucl_approx: approx nucleation temperature
                                                  nucl: nucleation temperature
                                                  perc: percolation temperature
                                                  compl: completion temperature
--epsturb=                              0.1       turbulence efficiency factor
                                                  >0: user defined
                                                  -1: upper bound
--pnlo_scaling=                         1         1 -> N NLO pressure
                                                  1: propto gamma
                                                  2: propto gamma^2
--checknlo=                             on        check for NLO stability
                                                  on: only keep NLO stable points
                                                  off: check disabled
--checkewsr=                            on        check for EWSR at high temperature
                                                  on: perform check and add info
                                                  keep_bfb: only keep BFB points
                                                  keep_ewsr: only keep EWSR points
                                                  off: check disabled
--maxpathintegrations=                  7         number of solutions of 1D equation =
                                                  number of path deformations + 1
--vevprofile=                           field     vev profile to solve BAU ('kink' or 'field')
--truncationscheme=                     minusvw   truncation scheme to be used
                                                  minusvw: R = -vw
                                                  <number>: R = number
                                                  variance: variance truncation
--moments=                              2         moments to solve the transport equations
--gwoutput=                             true      print GW wave output
--usegsl=                               true      use GSL library for minimization
--usecmaes=                             true      use CMAES library  for minimization
--usenlopt=                             true      use NLopt library for minimization
--usemultithreading=                    false     enable multi-threading for minimizers
--json=                                           use a json file instead of cli parameters

The following options for the Logger are available:

--logginglevel::default                 true      default output
--logginglevel::debug                   false     debug output
--logginglevel::disabled                          disable all output
--logginglevel::progdetailed            false     executable output
--logginglevel::minimizerdetailed       false     minimizer output
--logginglevel::transitiondetailed      false     transition calculation output
--logginglevel::mintracerdetailed       false     minimum tracking output
--logginglevel::bouncedetailed          false     bounce solution calculation output
--logginglevel::gwdetailed              false     gravitational wave calculation output
--logginglevel::vacuumprofiledetailed   false     vacuum profile calculation output
--logginglevel::fhckdetailed            false     baryogenesis using FHCK calculation output
--logginglevel::complete                false     enable all except minimizer output

The chosen Method for the thermal mass corrections is Arnold Espinosa
The implemented models are 
template
rxsm
cxsm
sm
n2hdm
cpinthedark
r2hdm
c2hdm

Not all required parameters are set.
