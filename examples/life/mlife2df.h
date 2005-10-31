C
      common /mlife2d/ opt_rows, opt_cols, opt_iter, opt_prows, 
     $                 opt_pcols
      common /mlife2dc/ opt_prefix
      integer opt_rows, opt_cols, opt_iter, opt_prows, opt_pcols
      character*80 opt_prefix
      integer BORN, DIES, MaxLRows, MaxLCols
      parameter (BORN=1)
      parameter (DIES=0)
      parameter (MaxLRows=100)
      parameter (MaxLCols=100)
