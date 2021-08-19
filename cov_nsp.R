


nsp1  = seq( 1, 180 )
nsp2  = seq( 181, 818 )
nsp3  = seq( 819, 2763 )
nsp4  = seq( 2764, 3263 ) 
nsp5  = seq( 3264, 3569 )
nsp6  = seq( 3570, 3859 )
nsp7  = seq( 3860, 3942 )
nsp8  = seq( 3943, 4140 )
nsp9  = seq( 4141, 4253 )
nsp10 = seq( 4254, 4392 )
nsp11 = seq( 4393, 4405 )
nsp12 = seq( 4393, 5324 )
nsp13 = seq( 5325, 5925 )
nsp14 = seq( 5926, 6452 )
nsp15 = seq( 6453, 6798 )
nsp16 = seq( 6799, 7096 )

name = paste0( "nsp", seq( 1,16 ) )

require( seqinr )

nsp = function( pos )
{
  start = 265
  S     = 21555
  
  if( pos < start | pos > S ){ stop( "Noncoding region" ) }
  
  pos_n = pos - start
  
  n_gene = 0
  while( pos_n > 0 )
  {
    n_gene   = n_gene + 1 # nsp?
    check_in = pos_n      # copy pos_n
    pos_n    = pos_n - length( get( name[n_gene] ) )*3
  }
  
  codon_ = check_in %% 3
  
  aa_n = ifelse( codon_ == 0, (check_in/3), (check_in %/% 3 + 1) )
  
  aa_name = paste0( "aa site - nsp ", n_gene )
  
  cat( paste0( aa_name, " ", aa_n, "\n" ) )
}



# nsp(3463)
# nsp(10132)



