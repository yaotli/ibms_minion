# simple script to annotate the N in final output 
# 
# package dependence: seqinr, stringr

suppressWarnings(suppressMessages( require( seqinr ) ) )
suppressWarnings(suppressMessages( require( stringr ) ) )


args = commandArgs( trailingOnly = TRUE )

if ( length(args) !=1 ) { stop( "try - Rscript post_artic.R folder", call. = FALSE ) } 

path_files = list.files( args[1], full.names = TRUE, recursive = TRUE )
sub_folder = gsub( "[^\\/]+$", "", grep( "\\.consensus\\.fasta", path_files, value = TRUE ) )

# message 1
mes = paste0( "...found ", length( sub_folder ), " samples\n" )
cat( mes )


for( i in 1: length( sub_folder ) )
{
  ### READIN ------
  
  path_sub_fd = list.files( sub_folder[i], full.names = TRUE )
  
  path_vcf1       = grep( "_1.vcf$", path_sub_fd, value = TRUE )
  path_vcf2       = grep( "_2.vcf$", path_sub_fd, value = TRUE )
  path_merged_vcf = grep( "merged.vcf$", path_sub_fd, value = TRUE )
  path_fail_vcf   = grep( "fail.vcf$", path_sub_fd, value = TRUE )
  path_primer_vcf = grep( "primers.vcf$", path_sub_fd, value = TRUE )
  path_cov_mask   = grep( "mask.txt$", path_sub_fd, value = TRUE )
  path_seq        = grep( "muscle.out.fasta$", path_sub_fd, value = TRUE )
  
  if( NA %in% c( path_vcf1, path_vcf2, path_merged_vcf, path_fail_vcf, path_primer_vcf, path_cov_mask, path_seq ) )
  {
    cat( paste0( "...incompleted outputs in subfolder ", sub_folder[i])  )
    next()
  }
  
  # vc1_2
  vcf1_vf2  = do.call( rbind, lapply( as.list( c( path_vcf1, path_vcf2 ) ), 
                                      function(x) tryCatch( read.table( file = x, 
                                                                        stringsAsFactors = FALSE, 
                                                                        colClasses = "character" ),
                                                            error = function(x) NULL ) ) )
  # fail
  fail_vcf = tryCatch( read.table( file = path_fail_vcf, stringsAsFactors = FALSE, colClasses = "character" ),
                       error = function(e){ return( NA ) } ) 
  
  if( is(fail_vcf, "data.frame") )
  {
    fail_vcf_2 = 
      do.call( rbind,
               lapply( as.list( fail_vcf$V8 ), 
                       function(x)
                       {
                         y = unlist( strsplit(x, ";") )
                         
                         z = str_match( y, "([^=]+)\\=([^=]+)" )
                         df = as.data.frame( t( as.matrix( z[,c(3)] ) ), stringsAsFactors = FALSE )
                         colnames( df ) = as.vector( t( as.matrix( z[,c(2)] ) ) )
                         
                         return( df )
                       } )  )
    
    fail_vcf = cbind( fail_vcf, fail_vcf_2 )
  }
  
  
  # merged
  merged_vcf = tryCatch( read.table( file = path_merged_vcf, stringsAsFactors = FALSE, colClasses = "character" ),
                         error = function(e){ return( NA ) } ) 
  
  if( is(merged_vcf, "data.frame") )
  {
    merged_vcf_2 = 
      do.call( rbind,
               lapply( as.list( merged_vcf$V8 ), 
                       function(x)
                       {
                         y = unlist( strsplit(x, ";") )
                         
                         z = str_match( y, "([^=]+)\\=([^=]+)" )
                         df = as.data.frame( t( as.matrix( z[,c(3)] ) ), stringsAsFactors = FALSE )
                         colnames( df ) = as.vector( t( as.matrix( z[,c(2)] ) ) )
                         
                         return( df )
                       } )  )
    
    merged_vcf = cbind( merged_vcf, merged_vcf_2 )
  }
  
  # primer
  primer_vcf = tryCatch(  read.table( file = path_primer_vcf, stringsAsFactors = FALSE, colClasses = "character" ),
                          error = function(e){ return( NA ) } ) 
  
  cov_mask  = read.table( path_cov_mask, stringsAsFactors = FALSE, header = FALSE )
  align_seq = seqinr::getSequence( seqinr::read.fasta( path_seq ) )

  
  ### CHECK CALLING METHODS ------
  
  .calling_method = 
    grepl( "nanopolish", read.delim( path_vcf1, sep = "\n", nrows = 2, header = FALSE )[2,1], ignore.case = TRUE )
  
  calling_method = ifelse( .calling_method, "nano", "medaka" )

    
  ### TERMINUS ------
  
  terminus    = unlist( apply( cov_mask[c(1, dim(cov_mask)[1]), -1], 1, function(x) seq( x[1], x[2] ) ) )
  OUTterminus = paste0( paste0( cov_mask[1, ]$V2, "-", cov_mask[1, ]$V3 ), "; ",
                        paste0( cov_mask[ dim(cov_mask)[1], ]$V2, "-", cov_mask[ dim(cov_mask)[1], ]$V3 ) )
  
  ### LOW-COVERAGE MASK ------
  
  n_lowCov_region = dim(cov_mask)[1] - 2
  
  #OUT_n_lowCov = sum( cov_mask[2:(n_lowCov_region+1), ]$V3 - cov_mask[2:(n_lowCov_region+1), ]$V2 ) + dim(cov_mask)[1]
  
  if( n_lowCov_region != 0 )
  {
    lowCov       = unlist( apply( cov_mask[2:(n_lowCov_region+1), -1], 1, function(x) seq( x[1], x[2] ) ) )
    OUT_lowCov   = paste0( apply( cov_mask[2:(n_lowCov_region+1), -1], 1, function(x) paste0( x, collapse = "-") ), 
                           collapse = "; " )
    
    Out_N_lowCov = length( lowCov )
    
  }else
  {
    lowCov       = NA
    OUT_lowCov   = ""
    Out_N_lowCov = 0
  }
  
  
  ### QUAL/TOTALREADS/SUPPORTFRACTION/HETEROGENOTYPE ------
  
  insert_i = which( align_seq[[2]] == "-")
  if( length( insert_i ) != 0 )
  {
    all_N = which( align_seq[[1]][ - insert_i ] == "n" )
    OUTinsertion = "True"
    
  }else { all_N = which( align_seq[[1]] == 'n' ); OUTinsertion = "False" }
  
  OUTall_N = length( all_N )
  
  inner_N    = setdiff( all_N, c( terminus, lowCov ) )
  OUTinner_N = length( inner_N )
  
  if( is(fail_vcf, "data.frame") )
  {
    n_vcf      = dim( fail_vcf )[1]
    n_vcf_note = rep( "", n_vcf ) 
    
    if( calling_method == "nano" )
    {
      for( v in 1:n_vcf ) 
      {
        supp_frac1 = as.numeric( strsplit( fail_vcf$SupportFractionByStrand[v], "," )[[1]][1] )
        supp_frac2 = as.numeric( strsplit( fail_vcf$SupportFractionByStrand[v], "," )[[1]][2] )
        
        if( as.numeric( fail_vcf[v, ]$TotalReads ) < 20 ){ n_vcf_note[v] = "Reads<20" 
        
        }else if ( ( as.numeric( merged_vcf[v, ]$V6 ) / as.numeric( fail_vcf[v, ]$TotalReads ) )  < 3 ){ n_vcf_note[v] = "low_QUAL" 
        
        }else if ( supp_frac1 < 0.5 | supp_frac2 < 0.5 ){ n_vcf_note[v] = "SuppFracStrand<0.5" 
        
        }else{ n_vcf_note[v] = "Inframe?" } 
      }
      
    }else
    {
      
      for( v in 1: n_vcf )
      {
        heter = gsub( "\\/", "", sapply( strsplit( fail_vcf$V10, ":" ), function(x) x[1] ) ) %in% c( "10", "01" )
        
        if( as.numeric( fail_vcf[v, ]$V6 )  < 20 ){ n_vcf_note[v] = "low_QUAL" 
        
        }else if( as.numeric( fail_vcf[v, ]$DP ) < 20 ){ n_vcf_note[v] = "Reads<20" 
        
        }else if( heter[v] ){ n_vcf_note[v] = "mixedGenotype" 
        
        }else{ n_vcf_note[v] = "Inframe?" } 
      }
    }
    
    
    ### UNKNOWN ------
    
    n_ref = nchar(fail_vcf$V4)
    n_alt = nchar(fail_vcf$V5)
    
    list_inner_N = as.list( as.numeric( fail_vcf$V2 ) )
    
    for( k in 1: length( n_ref ) )
    {
      if( n_ref[k]  >  n_alt[k]  )
      {
        list_inner_N[[ k ]] = seq( list_inner_N[[ k ]], ( list_inner_N[[ k ]] + n_ref[k]-1 ) )
        fail_vcf$V2[k]      = paste0( fail_vcf$V2[k], "Del", (n_ref[k]-n_alt[k]) )
      }
    }
    
    unknown_i = which( !inner_N %in% unlist( list_inner_N ) )
    
    if( length( unknown_i ) == 0 ){ OUTunknown = "" }else{ OUTunknown = paste0( inner_N[ unknown_i ], collapse = "; " ) }
    
  }else{ OUTunknown = "" }
  

  
  
  ### PRIMER REGIONS/LONELY PRESENT ------
  
  if( is( merged_vcf, "dataframe" ) )
  {
    not_merged_i = which( !vcf1_vf2$V2 %in% merged_vcf$V2 )
    
  }else{ not_merged_i = seq( 1, length( vcf1_vf2$V2 ) ) }
  
  if( length( not_merged_i ) == 0 ){ OUTnot_merged = "" }else{ OUTnot_merged = paste0( vcf1_vf2$V2[ not_merged_i ], collapse = "; " ) }
  
  if( is(primer_vcf, "data.frame") ){ OUTprimer = paste0( primer_vcf$V2, collapse = "; " ) }else{  OUTprimer = "" }
  
  
  
  ### SEQ DIFFERENCES ------
  
  mx = do.call( rbind, align_seq )
  
  var_site = 
    apply( mx, 2, 
           function(x)
           {
             y = unique(x)
             
             if( ( length( y ) > 1)  & ( !"n" %in% y  ) ) 
             {
               return( TRUE )
               
             }else{ FALSE }
           })
  
  OUTn_var = length( which( var_site ) )
  

  
  ### OUTPUT ------
  
  #1 touch the file 
  
  n_barcode = str_match( grep( "\\.consensus.fasta", path_sub_fd, value = TRUE ),
                         "barcode([0-9]+).consensus.fasta" )[,2]
  
  path_file = paste0( sub_folder[i], "barcode", n_barcode, ".report.txt" )
  
  if( TRUE %in% grepl( "barcode[0-9]+\\.report.txt", path_sub_fd ) )
  {
    cml_s = paste0( "echo '", "# ------------------------ START", "' > ", path_file )
    system( cml_s )
    
  }else
  {
    cml1 = paste0( "touch ", path_file )
    system( cml1 )  
    
    cml_s = paste0( "echo '", "# ------------------------ START", "' >> ", path_file )
    system( cml_s )
  }
  
  #2 print to the file
  
  l1 = paste0( "Total N positions: ", OUTall_N )
  l2 = paste0( "Total variants: ", OUTn_var )
  l3 = "# ------------------------"
  l4 = paste0( "Regions of low coverages: ", OUT_lowCov )
  l5 = paste0( "N in low coverages: ", Out_N_lowCov )    
  l6 = paste0( "Unsequenced terminus: ", OUTterminus )
  l7 = paste0( "N beyond terminus_lowCov: ", OUTinner_N )
  l8 = "# ------------------------"
  l9 = paste0( "Var not merged: ", OUTnot_merged )
  l10 = paste0( "Var at primer sites: ", OUTprimer )
  l11 = paste0( "Masked for unknown reason: ", OUTunknown )
  l12 = paste0( "Insersion: ", OUTinsertion )
  l13 = "# ------------------------"
  
  for(m in 1:13)
  {
    ll = get( paste0( "l", m ) )
    
    cml2 = paste0( "echo '", ll, "' >> ", path_file )
    system( cml2 )
  }
  
  if( is(fail_vcf, "data.frame") )
  {
    for(j in 1: length(n_vcf_note) )
    {
      cml3 = paste0( "echo '", fail_vcf$V2[j], ": ", n_vcf_note[j], "' >> ", path_file )
      system( cml3 )
    } 
  }else
  {
    cml3 = paste0( "echo '", "**Empty pass.vcf**", "' >> ", path_file )
    system( cml3 )  
  }
  
  cml_e = paste0( "echo '", "# ------------------------ COMPLETE", "' >> ", path_file )
  system( cml_e )
  
  # 
  cat("-")
}




#### VERSION ####
# 20210709

