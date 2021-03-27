Setting up parameters for simulation
====================================
Edit the *input.ini* or make your own *.ini* and run the code. The basic structure of *input.ini* is provided below,

.. code-block:: ini

  ;
  ; @file   input.ini
  ; @brief  PPDyn inputfile.
  ;
  scope = default

  [simbox]
  Lx  = 10.0    ; System length in X
  Ly  = 10.0    ; System length in Y
  Lz  = 10.0    ; System length in Z

  [particles]
  N     = 100     ; Number of particles
  Vxmax = 1.0     ; Maximum velocity in X
  Vymax = 1.0     ; Maximum velocity in Y
  Vzmax = 1.0     ; Maximum velocity in Z
  Temp  = 0.010   ;

  [boundary]
  btype = periodic ; Type of boundary

  [time]
  tmax  = 50.0    ; Final time
  dt    = 0.010   ; time step size

  [diagnostics]
  dumpPeriod  = 10    ; Data dump period
  dumpData    = False
  [options]
  parallelMode  = True  ;set to false to disable parallel
  ```
