# function to visualize variant positions
# 
# package dependence: seqinr

suppressWarnings(suppressMessages( require( seqinr ) ) )
suppressWarnings(suppressMessages( require( stringr ) ) )
suppressWarnings(suppressMessages( require( ggplot2 ) ) )
suppressWarnings(suppressMessages( require( ggrepel ) ) )
# install by: devtools::install_github("slowkow/ggrepel")


## FUNCTIONS ------

fastaEx = function(filedir = file.choose())
{
  file     <- seqinr::read.fasta(filedir, forceDNAtolower = FALSE)
  file_seq <- getSequence(file)
  file_seq <- lapply( file_seq, toupper )  
  file_id  <- attributes(file)$names
  
  return( list(seq = file_seq, 
               id  = file_id ) )
  #v201706
}

parseFEATURES = function( featTab = file.choose() )
{
  readin = read.delim( featTab, header = FALSE, fill = TRUE, skip = 1, stringsAsFactors = FALSE )
  
  n_row = nrow( readin )
  
  .class = c()
  .range = c()
  .name  = c()
  
  for( r in 1: n_row )
  {
    if( !readin$V3[r] == "" )
    {
      r_s = r
      
      if( !TRUE %in% (which( readin$V3 != "" ) > r) )
      {
        r_e = n_row
        
      }else
      {
        r_e = min( which( readin$V3 != "" )[ which( readin$V3 != "" ) > r ] ) - 1  
      }
      
      grid = readin[ r_s:r_e, ]
      
      .class = c( .class, grid$V3[ which( grid$V3 != "" ) ] )
      
      r_range = paste( grid$V1, grid$V2, sep = "-" )
      r_range = paste0( r_range[ !grepl( "^N", r_range) ], collapse = ";" )
      
      .range = c( .range, r_range )
      
       if( TRUE %in% grepl( "gene$|product$", grid$V4 ) )
       {
         r_name = grep( "gene$|product$", grid$V4 )
         .name  = c( .name, grid$V5[ r_name ]  )
         
       }else{ .name  = c( .name, ""  ) }
    }
  }
  
  df = data.frame( range = .range, class = .class, name = .name, stringsAsFactors = FALSE )
  
  return( df )
}

linkSingle = function( v )
{
  if( length(v) > 0 )
  {
    v = sort(v)
    
    v_point = c(1,  ( which( c( v[-1] - v[-length(v)] ) != 1 ) + 1 ) )
    
    out = list()
    for( j in 1: length( v_point ) )
    {
      if( j == length( v_point ) )
      {
        out[[j]] = v[ v_point[j]: length(v) ]
      }else
      {
        out[[j]] = v[ v_point[j]: (v_point[j+1]-1) ]  
      }
    } 
    return( out )
    
  }else
  {
    return( NA )  
  }
}

sendLines = function( lines, type = "v" )
{
  if( type == "v")
  {
    for( l in 1: length( lines ) )
    {
      cat( paste0( lines[l], " \n" ) )
    }
  }else
  {
    cat( lines )
    cat( "\n" )
  }
}

## VARVIZ ------

varViz = function( infile        = file.choose(), 
                   batchName     = "",
                   refSeq        = 1,             # the first seq in .fas is NC_045512 or MN908947
                   showType      = "A",           # show [A]mino or [N]ucleotide
                   region        = "S",           # coding region; only if [A] mode is used 
                   saveFile      = FALSE,         # auto save a table
                   savePlot      = FALSE,         # auto save a pdf figure 
                   #ignoreN      = TRUE,          # FIGURE: not showing N 
                   ignoreTermini = 70 )           # FIGURE: length of unsequenced nucleotide termini ignored
{
  # NC_045512 = parseFEATURES()
  # NC_045512 = NC_045512[ which( NC_045512$class == "CDS" ), ]
  NC_045512.range =
    c(
      "266-13468;13468-21555", "266-13483",
      "21563-25384",           "25393-26220",          
      "26245-26472",           "26523-27191",          
      "27202-27387",           "27394-27759",          
      "27756-27887",           "27894-28259",          
      "28274-29533",           "29558-29674" )
  NC_045512.name = 
    c(
      "ORF1ab polyprotein",          "ORF1a polyprotein",          
      "surface glycoprotein",        "ORF3a protein",              
      "envelope protein",            "membrane glycoprotein",      
      "ORF6 protein",                "ORF7a protein",              
      "ORF7b",                       "ORF8 protein",               
      "nucleocapsid phosphoprotein", "ORF10 protein" )
  
  
  region_range = 
    sapply( strsplit( NC_045512.range, ";"  ), 
            function(x)
            {
              y = c()
              
              for( j in 1: length(x) )
              {
                s = as.numeric( str_match( x[j], "([0-9]+)-([0-9]+)" )[,2] )
                e = as.numeric( str_match( x[j], "([0-9]+)-([0-9]+)" )[,3] )
                
                y = c( y, seq(s,e) ) 
              }
              return( y )
            })
  
  region_names = c( "ORF1ab", "ORF1a", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10" )
  
  
  ## INPUT ------
  
  seq = fastaEx( infile )$seq
  id  = fastaEx( infile )$id
  
  print( id )
  
  ## * readline 1 ------
  
  sub1 = readline( paste0( "Choose 1st seq [1-", length(id), "], enter [q] to quit\n" ) )
  
  while( !( (sub1 %in% as.character( seq( 1: length(id) ) ) ) | (sub1 == "q") ) )
  {
    sub1 = readline( paste0( "Choose 1st seq [1-", length(id), "], enter [q] to quit\n" ) )
  }
  
  if( sub1 == "q" ){ stop( "Quit\n" ) }
  
  ## * readline 2 ------
  
  sub2 = readline( paste0( "Choose 2nd seq [1-", length(id), "], enter [q] to quit\n" ) )
  
  while( !( (sub2 %in% as.character( seq( 1: length(id) ) ) ) | (sub2 == "q") ) | (sub1 == sub2) )
  {
    sub2 = readline( paste0( "Choose 2nd seq [1-", length(id), "], enter [q] to quit\n" ) )
  }
  
  if( sub2 == "q" ){ stop( "Quit\n" ) } 
  
  
  ## N VARIANT ------ 
  
  # 
  n.pos  = c()
  n.char = c()
  n.ref  = c()
  n.sub1 = c()
  n.sub2 = c()
  #
  a.pos  = c()
  a.char = c()
  a.ref  = c()
  a.sub1 = c()
  a.sub2 = c()
  #
  
  
  mx  = do.call( rbind, seq[ c( refSeq, as.numeric( c( sub1, sub2 ) ) ) ] )
  lth = dim(mx)[2]
  
  ## * Insersion ------
  
  insert = which(mx[1,]== "-")
  if( length( insert ) > 0 )
  {
    sum_gap   = 0 
    ls_insert = linkSingle( insert )
  
    for( j in 1: length( ls_insert ) )
    {
      in_s = ls_insert[[j]][1] -1
      in_e = ls_insert[[j]][ length( ls_insert[[j]] ) ]
      
      n.pos  = c( n.pos, (in_s- sum_gap) )
      n.char = c( n.char, paste0("I", (in_e - in_s) ) )
      n.ref  = c( n.ref, c2s( mx[1, in_s:in_e] ) )
      n.sub1 = c( n.sub1, c2s( mx[2, in_s:in_e] ) )
      n.sub2 = c( n.sub2, c2s( mx[3, in_s:in_e] ) )
      
      sum_gap = sum_gap + length( ls_insert[[j]] ) 
    }
    
    mx = mx[, -insert]
  }
  
  if( !dim(mx)[2] == 29903 ){ stop( "Insertion" ) }
  
  ## * Deletion ------
  
  list_del = lapply( list( mx[ 2, ], mx[ 3, ]), function(x) which(x == "-" ) )
  
  if( ( length( list_del[[1]] != 0 ) ) | ( length( list_del[[2]] != 0 ) ) )
  {
    ls_del = lapply( list_del, linkSingle )
    
    likely_na = sapply( ls_del, function(x){ (!is.list(x) ) & (length(x)==1) } )
    
    if( TRUE %in% likely_na ){ ls_del = ls_del[ - which( likely_na ) ] }
    
    for( j in 1: length( ls_del ) )
    {
      for( k in 1: length( ls_del[[j]] ) )
      {
        del_s = ls_del[[j]][[k]][1]
        del_e = ls_del[[j]][[k]][ length( ls_del[[j]][[k]] ) ]
        
        n.pos  = c( n.pos, del_s )
        n.char = c( n.char, paste0("D", (del_e - del_s +1 ) ) )
        n.ref  = c( n.ref, c2s( mx[1, del_s:del_e] ) )
        n.sub1 = c( n.sub1, c2s( mx[2, del_s:del_e] ) )
        n.sub2 = c( n.sub2, c2s( mx[3, del_s:del_e] ) ) 
      }
    }
    
  }
  
  
  ## * N ------
  
  list_N = lapply( list( mx[ 2, ], mx[ 3, ]), function(x) which(x == "N" ) )
    
  if( ( length( list_N[[1]] != 0 ) ) | ( length( list_N[[2]] != 0 ) ) )
  {
    for( j in 1: 2 )
    {
      if( length( list_N[[j]] ) == 0 ){ next() }
      
      n.pos  = c( n.pos, list_N[[j]] )
      n.char = c( n.char, rep( "N", length( list_N[[j]] ) ) )
      n.ref  = c( n.ref, mx[1, list_N[[j]] ] )
      n.sub1 = c( n.sub1, mx[2, list_N[[j]] ] )
      n.sub2 = c( n.sub2, mx[3, list_N[[j]] ] )  
    }
  }
  
  
  ## * Mutation ------
  
  var_site = 
    apply( mx, 2, 
           function(x)
           {
             y = unique(x)
             
             if(  (length( y ) > 1) & ( !"N" %in% x ) & ( !"-" %in% x ) )
             {
               return( TRUE )
               
             }else{ FALSE }
           })
  
  var_site_i = which( var_site )
  
  n.pos  = c( n.pos, var_site_i )
  n.char = c( n.char, rep( "M", length(var_site_i) ) )
  n.ref  = c( n.ref, mx[1, var_site_i ] )
  n.sub1 = c( n.sub1, mx[2, var_site_i ] )
  n.sub2 = c( n.sub2, mx[3, var_site_i ] )  
  
  ## COMBINE ------
  
  nu_df = data.frame( pos = n.pos, char = n.char, ref = n.ref, sub1 = n.sub1, sub2 = n.sub2, stringsAsFactors = F )
  
  ## AMINO ACID ------
  
  gene_i = match( region, region_names )
  
  if( (showType == "A") & ( !is.na( gene_i ) ) )
  {
    range_i = region_range[[ gene_i ]]
    
    mx_sub = mx[ ,range_i ]
    if( (dim(mx_sub)[2] %%3) != 0 ){ warning( "...may not inframe \n" ) }
    
    mx_aa = t( apply( mx_sub, 1, translate ) )
    
    ## * Del_Outframe ------
    
    list_x = lapply( list( mx_aa[ 2, ], mx_aa[ 3, ]), function(x) which(x == "X" ) )
    
    if( ( length( list_x[[1]] != 0 ) ) | ( length( list_x[[2]] != 0 ) ) )
    {
      ls_x = lapply( list_x, linkSingle )
      
      likely_na = sapply( ls_x, function(x){ (!is.list(x) ) & (length(x)==1) } )
      
      if( TRUE %in% likely_na ){ ls_x = ls_x[ - which( likely_na ) ] }
      
      for( j in 1: length( ls_x ) )
      {
        for( k in 1: length( ls_x[[j]] ) )
        {
          x_s = ls_x[[j]][[k]][1]
          x_e = ls_x[[j]][[k]][ length( ls_x[[j]][[k]] ) ]
          
          a.pos  = c( a.pos, x_s )
          a.char = c( a.char, paste0("X", (x_e - x_s + 1) ) )
          a.ref  = c( a.ref, c2s( mx_aa[1, x_s:x_e] ) )
          a.sub1 = c( a.sub1, c2s( mx_aa[2, x_s:x_e] ) )
          a.sub2 = c( a.sub2, c2s( mx_aa[3, x_s:x_e] ) ) 
        }
      }
      
    }
    
    ## * Variant ------
    
    var_site = 
      apply( mx_aa, 2, 
             function(x)
             {
               y = unique(x)
               
               if(  (length( y ) > 1) & ( !"X" %in% y ) )
               {
                 return( TRUE )
                 
               }else{ FALSE }
             } )
    
    var_site_i = which( var_site )
    
    a.pos  = c( a.pos, var_site_i )
    a.char = c( a.char, rep("V", length( var_site_i ) ) )
    a.ref  = c( a.ref, mx_aa[1, var_site_i] ) 
    a.sub1 = c( a.sub1, mx_aa[2, var_site_i] ) 
    a.sub2 = c( a.sub2, mx_aa[3, var_site_i] )  
    
    aa_df = data.frame( pos = a.pos, char = a.char, ref = a.ref, sub1 = a.sub1, sub2 = a.sub2, stringsAsFactors = F )
    
    ## REPORT AA ------
    
    # V ---
    v_i = intersect( which( ! aa_df$sub1 == aa_df$sub2 ), which( aa_df$char == "V" ) )
    
    if( length( v_i ) > 0 )
    {
      out_v     = paste0( paste0( region, ": " ), aa_df$sub1[ v_i ], " ", aa_df$pos[ v_i ], " ", aa_df$sub2[ v_i ]  )
      out_v_pos = aa_df$pos[ v_i ] 
      
    }else
    {
      out_v = ""
    }
    
    # X ---
    x_i = intersect( which( ! aa_df$sub1 == aa_df$sub2 ), grep( "^X", aa_df$char )  )
    
    if( length( x_i ) > 0 )
    {
      out_x     = paste0( paste0( region, ": " ), aa_df$pos[ x_i ], " ", aa_df$char[ x_i ] )
      out_x_pos = aa_df$pos[ x_i ]
      
    }else
    {
      out_x = ""
    }
  }
  
  
  ## REPORT NT ------
  
  # N ---
  out_num_N1 = paste0( length( which( nu_df$sub1 == "N" ) ), " N in seq 1" )
  out_num_N2 = paste0( length( which( nu_df$sub2 == "N" ) ), " N in seq 2" )

  # M ---
  
  m_i = intersect( which( ! nu_df$sub1 == nu_df$sub2 ), which( nu_df$char == "M" ) )
  
  if( length( m_i ) > 0 )
  {
    out_m     = paste0( nu_df$sub1[ m_i ], " ", nu_df$pos[ m_i ], " ", nu_df$sub2[ m_i ]  )
    out_m_pos = nu_df$pos[ m_i ]
    
  }else
  {
    out_m = ""
  }
  
  # INDEL ---
  indel_i = intersect( which( ! nu_df$sub1 == nu_df$sub2 ), which( grepl( "^D|^I", nu_df$char ) ) )
  if( length( indel_i ) > 0 )
  {
    out_ind = paste0( nu_df$pos[ indel_i ], " ", nu_df$char[ indel_i ]  )
    out_ind = gsub( "I", "In", out_ind )
    out_ind = gsub( "D", "Del", out_ind )
    
    out_ind_pos = nu_df$pos[ indel_i ]
    
  }else
  {
    out_ind = ""
  }
  
  
  ## OUT ------
  
  cat( "------ note ------------ \n" )
  sendLines( c( out_num_N1, out_num_N2 ) )
  
  if( exists( "out_v_pos" ) )
  {
    cat( "------ Protein ------------ \n" )
    
    sendLines( out_v[ order(out_v_pos) ] )
  }
  
  if( exists( "out_x_pos" ) )
  {
    cat( "------ [Pos] [X][range] ------------ \n" )
    sendLines( out_x[ order(out_x_pos) ] )
  }
  
  if( exists( "out_m_pos" ) )
  {
    cat( "------ Nucleotide ------------ \n" )
    
    sendLines( out_m[ order(out_m_pos) ], "a" )
  }
  
  if( exists( "out_ind_pos" ) )
  {
    cat( "------ [Pos] [Del/In][range] ------------ \n" )
    sendLines( out_ind[ order(out_ind_pos) ] )
  }

  ## SAVE FILE ------
  
  if( saveFile & (dim(nu_df)[1] != 0) )
  {
    #nu_df = nu_df[ -which( nu_df$char == "N" ), ]
      
    if( !exists( "aa_df" ) )
    {
      out_df = data.frame( nu_df, type = rep( "nt", dim(nu_df)[1] ), stringsAsFactors = FALSE )
      
    }else
    {
      out_df =   
        rbind( data.frame( nu_df, type = rep( "nt", dim(nu_df)[1] ), stringsAsFactors = FALSE ), 
               data.frame( aa_df, type = rep( paste0( "aa: ", region ), dim(aa_df)[1] ), stringsAsFactors = FALSE ) )
    }
    
    att_df = out_df[1,]
    att_df[,1] = 0
    att_df[,2] = paste0( "FILE= ", infile, "; SEQ1= ", id[ as.numeric(sub1) ], "; SEQ2= ", id[ as.numeric(sub2) ] )
    att_df[,3:6] = ""
    
    out_df = rbind( out_df, att_df )
    
    filename = gsub( ".fasta", paste0( ".", Sys.time(), ".tsv"), infile )
    
    write.table( out_df, sep = "\t", file = filename, quote = FALSE, row.names = FALSE )
  }
  
  
  ## VISUAL PREP ------
  
  # env ---
  
  lth_e = ifelse( showType == "A", length( range_i)/3+1, 29903 )
  
  df_fig_line = data.frame( h = c(1, 1), lth = c(1, lth_e ), seg = c( "all", "all" ) )
  
  if( (region == "S") & ( showType == "A") )
  {
    df_fig_line = rbind( df_fig_line, 
                         data.frame( h = 1, lth = c( 13, 685, 686, lth_e ), seg = c( "S1", "S1", "S2", "S2" ) ) )
  }
  
  # n ---
  if( length( which( nu_df$char == "N" ) ) > 0 )
  {
    n_pos = data.frame( pos = nu_df$pos[ which( nu_df$char == "N" ) ], h = 1, anno = "N", type = "n" )
    
  }else
  {
    n_pos = data.frame( pos = 0, h = 1, anno = "N", type = "n" )  
  }
  
  # m ---
  if( length( m_i ) > 0 )
  {
    m_pos = data.frame( pos  = nu_df$pos[ m_i ], h = 1, 
                        anno = paste0( nu_df$sub1[ m_i ], nu_df$pos[ m_i ], nu_df$sub2[ m_i ] ),
                        type = "m")
  }else
  {
    m_pos = data.frame( pos = 0, h = 1, anno = "NA", type = "m" )    
  }
  
  # indel ---
  if( length( indel_i ) > 0 )
  {
    indel_pos = data.frame( pos = nu_df$pos[ indel_i ], h = 1, 
                            sub1 = nu_df$sub1[ indel_i ], sub2 = nu_df$sub2[ indel_i ], anno = nu_df$char[ indel_i ],
                            type = "indel",
                            stringsAsFactors = FALSE )
    
    for( x in 1: length( indel_pos$anno ) )
    {
      if( grepl( "^I", indel_pos$anno[x] ) )
      {
        target            = ifelse( grepl( "-$", indel_pos$sub1[x] ), 2, 1 )
      }else
      {
        target            = ifelse( grepl( "^-", indel_pos$sub1[x] ), 1, 2 )
      } 
      
      .char = str_match( indel_pos$anno[x], "(^[DI])([0-9]+)" )[,2]
      .lth  = as.numeric( str_match( indel_pos$anno[x], "(^[DI])([0-9]+)" )[,3] )
      
      indel_pos$anno[x] = paste0( #"S", target, ": ", 
                                  .char, indel_pos$pos[x], "-", (indel_pos$pos[x]+.lth-1)  )
    }
    
    indel_pos = indel_pos[ -c(3,4) ]
    
  }else
  {
    indel_pos = data.frame( pos = 0, h = 1, anno = "NA", type = "indel" )    
  }
  
  # v & X ---
  if( (showType == "A") & ( !is.na( gene_i ) ) )
  {
    # v ---
    if( length( v_i ) > 0 )
    {
      v_pos = data.frame( pos  = aa_df$pos[ v_i ], h = 1, 
                          anno = paste0( aa_df$sub1[ v_i ], aa_df$pos[ v_i ], aa_df$sub2[ v_i ] ),
                          type =  "v" )
    }else
    {
      v_pos = data.frame( pos = 0, h = 1, anno = NA, type = "v" )  
    }
    
    # x ---
    if( length( x_i ) > 0 )
    {
      x_pos = data.frame( pos  = aa_df$pos[ x_i ], h = 1, 
                          sub1 = aa_df$sub1[ x_i ], sub2 = aa_df$sub2[ x_i ], anno = aa_df$char[ x_i ],
                          type = "X", 
                          stringsAsFactors = FALSE )
      
      for( x in 1: length( x_pos$anno ) )
      {
        target = ifelse( grepl( "^X", x_pos$sub1[x] ), 1, 2 )
        
        .char = str_match( x_pos$anno[x], "(^[X])([0-9]+)" )[,2]
        .lth  = as.numeric( str_match( x_pos$anno[x], "(^[X])([0-9]+)" )[,3] )
        
        x_pos$anno[x] = paste0( #"S", target, ": ", 
                                .char, x_pos$pos[x], "-", (x_pos$pos[x]+.lth-1)  )
      }
      
      x_pos = x_pos[ -c(3,4) ]
      
    }else
    {
      x_pos = data.frame( pos = 0, h = 1, anno = "NA", type = "X" )    
    }
  }
  
  # combine ---
  if( (showType == "A") & ( !is.na( gene_i ) ) )
  {
    nt_com_pos = rbind( n_pos, m_pos, indel_pos )
    nt_com_pos = nt_com_pos[ which( nt_com_pos$pos != 0 ), ]
    
    aa_com_pos = rbind( v_pos, x_pos )
    aa_com_pos = aa_com_pos[ which( aa_com_pos$pos != 0 ), ]
    
  }else
  {
    nt_com_pos = rbind( n_pos, m_pos, indel_pos )
    nt_com_pos = nt_com_pos[ which( nt_com_pos$pos != 0 ), ]  
  }
  
  ## * IGNORE N ------
  if( ignoreTermini != 0 )
  {
    ls_ignoreTermini = c( seq( 1, ignoreTermini ), seq( (29903-ignoreTermini+1), 29903 ) )
    
    nt_com_pos = nt_com_pos[ !nt_com_pos$pos %in% ls_ignoreTermini, ]
  }
  
  ## VISUAL ------
  
  if( showType == "A" )
  {
    fig = 
      ggplot( ) + 
      geom_line( data =  subset( df_fig_line, seg == "S1" ) , aes( x = lth, y = h, group = seg ), size = 4, color = "#d62728", alpha = 0.5 ) +
      geom_line( data =  subset( df_fig_line, seg == "S2" ) , aes( x = lth, y = h, group = seg ), size = 4, color = "#1f77b4", alpha = 0.5 ) +
      geom_line( data =  subset( df_fig_line, seg == "all" ) , aes( x = lth, y = h, group = seg ), size = 1.5 ) +
      geom_point( data = subset( aa_com_pos, type == "v" ), aes( x = pos, y = h ), color = "#d62728", size = 0.5 ) +
      geom_text_repel( data = subset( aa_com_pos, type == "v" ), aes( x = pos, y = h, label = anno ), 
                       size = 2.5,
                       force             = 0.5,
                       nudge_y           = 0.025,
                       direction         = "x",
                       hjust             = 0,
                       segment.size      = 0.2,
                       angle             = 90) + 
      
      geom_text_repel( data = subset( aa_com_pos, type == "X" ), aes( x = pos, y = h, label = anno ), 
                       size = 2.5,
                       force             = 0.5,
                       nudge_y           = -0.01,
                       direction         = "x",
                       vjust             = 0,
                       segment.size      = 0.2,
                       angle             = 90 ) + 
      
      scale_y_continuous( breaks = 1, label = ifelse( batchName == "" , region, batchName ),
                          limits = c(0.95, 1.05) ) +
      scale_x_continuous( breaks = c( seq( 0, lth_e, by = 500 ), lth_e ), 
                          label  = c( seq( 0, lth_e, by = 500 ), lth_e ) ) +
      theme_bw() +
      theme( panel.grid.minor   = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.border       = element_blank(),
             axis.title         = element_blank(),
             axis.ticks         = element_blank() )
    
  }else
  {
    fig = 
      ggplot( ) + 
      geom_line( data = df_fig_line, aes( x = lth, y = h, group = seg ), size = 1.5 ) + 
      geom_point( data = subset( nt_com_pos, type == "n" ), aes( x = pos, y = h ),  color = "#7f7f7f", shape = "I", size = 2 ) +
      geom_point( data = subset( nt_com_pos, type == "m" ), aes( x = pos, y = h ), color = "#d62728", size = 0.5 ) +
      geom_text_repel( data = subset( nt_com_pos, type == "m" ), aes( x = pos, y = h, label = anno ), 
                       size = 2,
                       force             = 0.5,
                       nudge_y           = 0.025,
                       direction         = "x",
                       hjust             = 0,
                       segment.size      = 0.2,
                       angle             = 90) + 
      
      geom_text_repel( data = subset( nt_com_pos, type == "indel" ), aes( x = pos, y = h, label = anno ), 
                       size = 2,
                       force             = 0.5,
                       nudge_y           = -0.01,
                       direction         = "x",
                       vjust             = 0,
                       segment.size      = 0.2,
                       angle             = 90 ) + 
      
      scale_y_continuous( breaks = 1, label = ifelse( batchName == "" , "Gene", batchName ),
                          limits = c(0.95, 1.05) ) +
      theme_bw() +
      theme( panel.grid.minor   = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.border       = element_blank(),
             axis.title         = element_blank(),
             axis.ticks         = element_blank() )  
    
    
    
    
  }
  
  print( fig )
  
  if( savePlot )
  {
    filename = gsub( ".fasta", paste0( ".", Sys.time(), ".pdf"), infile )
    
    suppressWarnings( ggsave( plot = fig, filename = filename, width = 5, height = 2, units = "in" ) )
  }
  
}



## RUN ------

varViz( region = "S", savePlot = T, ignoreTermini = 65 )






#### VERSION ####
# 20210721



