#' Construct a lolliplot
#'
#' Given a data frame construct a plot displaying mutations on a transcript
#' framework.
#' @name lolliplot
#' @param x Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name",
#' "gene", and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p.notation (see details).
#' @param y Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name"
#' and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p. notation (optional, see details).
#' @param z Object of class data frame with rows representing regions of
#' interest. The data frame must contain columns with the following names
#' "description", "start", "stop" (optional see details).
#' @param fillCol Character string specifying the column name of the argument
#' supplied to parameter x on which to colour the lollis representing mutations
#' (see details).
#' @param labelCol Character string specifying the column name of the argument
#' supplied to parameter x from which to extract and display text corresponding
#' to mutations (see details).
#' @param txtAngle Integer specifying the angle of label text to be plotted if
#' an argument is supplied to the labelCol parameter.
#' @param txtSize Integer specifying the size of label text to be plotted if an
#' argument is supplied to the labelCol parameter.
#' @param pntSize Integer specifying the size of lolli points representing
#' mutations.
#' @param proteinColour Character string specifying the background colour of the
#' protein.
#' @param obsA.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter x (see details and vignette).
#' @param obsA.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter x (see details and vignette).
#' @param obsB.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter y (see details and vignette).
#' @param obsB.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter y (see details and vignette).
#' @param sideChain Boolean specifying if amino acid sidechain data should be
#' plotted in lieu of protein domains (see details).
#' @param species A valid species from which to retrieve protein domain and
#' sequence data for a given transcript (see details).
#' @param maxLolliStack Integer specifying the cutoff for the maximum number of
#' lollis allowed to be stacked at a single position.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param paletteA Character vector specifying colours for protein domains,
#' valid only if sideChain==FALSE.
#' @param paletteB Character vector specifying colours for lollis representing
#' mutations, valid only if argument is supplied to fillCol.
#' @param host Host to connect to for biomaRt queries (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details lolliplot is a function designed to display mutation information in
#' the context of a protien identified by an ensembl transcript id. The
#' lolliplot function will query ensembl via biomart to retrieve sequence and
#' domain information in order to construct a representation of a protein and
#' therefore requires an internet connection. A value must be supplied to the
#' species parameter (defaults to hsapiens) in order for a successful biomart
#' query. Valid arguments to this field are those species with datasts available
#' via ensembl. please specify species in lowercase without a period
#' (i.e. hsapiens instead of H.sapiens), lolliplot will inform the user of
#' available species if input to the species parameter is not recognized.
#' Further lolliplot will build a protein framework based on sequence data
#' obtained from biomaRt, by default this will default to the latest ensembl
#' version. In order for the most accurate representation the annotation version
#' of the mutations given to lolliplot should match the annotation version used 
#' by biomaRt. The annotation version used by biomaRt can be changed via the 
#' host paramter (see vignette for more details).
#' 
#' lolliplot is capable of plotting two seperate sets of data on the protein
#' representation specified by parameters `x` and `y`, the data supplied to
#' these parameters will be plotted on the top and bottom of the protein
#' respectively. Note that input to these parameters is expected to correspond
#' to a single ensembl transcript and that values in the "amino_acid_change"
#' columns are required to be in p. notation (i.e. p.V600E). Further lolliplot
#' is able to plot custom domain annotation if supplied via the parameter `z`,
#' this will override domain information obtained from biomart.
#' 
#' lolliplot uses a forcefield model from the package FField to attract and 
#' repulse lollis. The parameters for this force field model are set to
#' reasonable defaults however may be adjusted via the obsA... and obsB...
#' family of parameters. Please see the package FField available on cran for
#' a description of these parameters. Note that the time to construct the
#' lolliplot will in large part depend on the number of mutations and the values
#' supplied to the forcefield parameters.
#' 
#' @examples
#' # Create input data
#' data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
#' data <- as.data.frame(cbind(data, 'ENST00000269305'))
#' colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')
#'
#' # Call lolliplot
#' lolliplot(data)
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @export

lolliplot <- function(x, y=NULL, z=NULL, fillCol=NULL, labelCol=NULL,
                      txtAngle=45, txtSize=5, pntSize=4,
                      proteinColour='#999999', obsA.rep.fact=5000,
                      obsA.rep.dist.lmt=500, obsA.attr.fact=.1, obsA.adj.max=.1,
                      obsA.adj.lmt=.5, obsA.iter.max=50000, obsB.rep.fact=5000,
                      obsB.rep.dist.lmt=500, obsB.attr.fact=.1, obsB.adj.max=.1,
                      obsB.adj.lmt=.5, obsB.iter.max=50000,
                      sideChain=FALSE, species="hsapiens",
                      maxLolliStack=NULL, plotLayer=NULL, paletteA=NULL,
                      paletteB=NULL, host="www.ensembl.org", out="plot")
{
  # Perform quality check
  input <- lolliplot_qual(x, y, z)
  x <- input[[1]]
  y <- input[[2]]
  z <- input[[3]]
  
  # extract transcript id and subset data y on that id if it exists
  transcriptID <- as.character(x$transcript_name[1])
  if(!is.null(y))
  {
    y <- y[y$transcript_name == transcriptID,]
  }
  
  # extract HUGO gene name
  gene <- as.character(x$gene[1])
  
  # Obtain length of protein
  result <- lolliplot_transcriptID2codingSeq(transcriptID,
                                             species=species,
                                             host=host)
  codingSeq <- result$coding
  cdsLen <- result$cds_length
  
  
  # Get the sequence length in AA, perform quality checks along the way
  residueSeq <- lolliplot_DNAconv(codingSeq, to="residue")    
  # If it is requested grab the sidechain information and bind to residues
  if(sideChain==TRUE)
  {
    sidechain <- lolliplot_DNAconv(codingSeq, to="sidechain")
    AAsequence <- cbind(sidechain, residueSeq)
    AAsequence <- as.data.frame(AAsequence)
    AAsequence$coord <- seq(from=1, to=nrow(AAsequence))
  } else {
    AAsequence <- NULL
  }
  
  # if there are any stop codons remove them as they are not considered part
  # of the protein
  if(any(residueSeq %in% c("OPAL", "OCHRE", "AMBER")))
  {
    stopRes <- c("OPAL", "OCHRE", "AMBER")
    residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
    if(!is.null(AAsequence))
    {
      AAsequence <- AAsequence[-which(AAsequence$residueSeq %in% stopRes),]
    }
  }
  
  # grab the length of the protein in Amino Acids
  proteinLength <- length(residueSeq)    
  
  # if z is specified plot that instead of fetching the domain information
  if(!is.null(z))
  {
    geneData <- lolliplot_constructGene(gene, z, proteinLength)
  } else {
    # extract protien domain data
    protein_domain <- lolliplot_fetchDomain(transcriptID,
                                            species=species,
                                            host=host)
    
    # construct gene from data collected
    geneData <- lolliplot_constructGene(gene, protein_domain, proteinLength)
  }
  
  # construct data frame of observed mutations for top track
  observed_mutation <- lolliplot_mutationObs(x, 'top', fillCol, labelCol,
                                             obsA.rep.fact, obsA.rep.dist.lmt,
                                             obsA.attr.fact, obsA.adj.max,
                                             obsA.adj.lmt, obsA.iter.max)
  observed_mutation <- lolliplot_reduceLolli(observed_mutation,
                                             maxLolliStack)
  
  # construct data frame of observed mutations for bottom track
  if(!is.null(y))
  {
    observed_mutation2 <- lolliplot_mutationObs(y, 'bottom', fillCol,
                                                labelCol, obsB.rep.fact,
                                                obsB.rep.dist.lmt,
                                                obsB.attr.fact,
                                                obsB.adj.max, obsB.adj.lmt,
                                                obsB.iter.max)
    observed_mutation2 <- lolliplot_reduceLolli(observed_mutation2,
                                                maxLolliStack)
  } else {
    observed_mutation2 <- NULL
  }
  
  # construct the lolliplot
  plot <- lolliplot_buildMain(geneData, length, observed_mutation,
                              observed_mutation2,fillCol, labelCol,
                              txtAngle, txtSize, pntSize,
                              proteinColour, AAsequence,
                              plot_sidechain=sideChain, layers=plotLayer,
                              paletteA=paletteA, paletteB=paletteB)
  
  # Decide what to output
  dataOut <- list("gene"=geneData,
                  "observed_mutation"=observed_mutation,
                  "observed_mutation2"=observed_mutation2)
  output <- multi_selectOut(data=dataOut, plot=plot, out=out)
  return(output)
}

#' Convert AA to side chain classification
#' 
#' Given the 1 letter code an amino acid, return the side chian classification
#' @name lolliplot_AA2sidechain
#' @param x Character of length 1 giving the 1 letter amino acid code
#' @return Object of class character

lolliplot_AA2sidechain <- function(x)
{
  # Coerce all AA changes to uppercase and then apply switch statement
  x <- toupper(x)
  x <- switch(EXPR=x, "F"="Nonpolar", "L"="Nonpolar", "S"="Polar",
              "Y"="Polar", "C"="Polar", "W"="Nonpolar", "L"="Nonpolar",
              "P"="Nonpolar", "H"="Basic", "Q"="Polar", "R"="Basic",
              "I"="Nonpolar", "M"="Nonpolar", "T"="Polar", "N"="Polar",
              "K"="Basic", "S"="Polar", "R"= "Basic", "V"="Nonpolar",
              "A"="Nonpolar", "D"= "Acidic", "E"="Acidic", "G"="Polar")
  
  return(x)
}

#' Construct Lolliplot
#'
#' Construct Lolliplot given gene and mutation data
#' @name lolliplot_buildMain
#' @param gene_data object of class dataframe giving protien domain and gene
#' information
#' @param length integer specifying the length of the protien in amino acids
#' @param mutation_observed object of class data frame specifying mutations
#' observed in input file
#' @param mutation_observed2 optional object of class data frame specifying
#' additional mutations for bottom track
#' @param fill_value character string specifying the column on which to colour
#' mutation points
#' @param label_column character string specifying the column containing the
#' labels to attach to mutation points
#' @param plot_text_angle numeric value specifying the angle of text to be
#' plotted
#' @param plot_text_size numeric value specifying the size of text to be plotted
#' @param point_size numeric value specigying the size of mutation points
#' @param gene_colour color to shade plotted gene
#' @param sequence_data object of class dataframe giving AA sequence, sidechain,
#' and coord required if plot_sidechain is true
#' @param plot_sidechain boolean specifying whether to plot the AA sidechain
#' instead of domain information
#' @param layers additional ggplot2 layers to plot
#' @param paletteA Character vector specifying colours for gene features
#' @param paletteB Character vector specifying colours for lolli features
#' @return a ggplot2 object
#' @import ggplot2

lolliplot_buildMain <- function(gene_data, length, mutation_observed,
                                mutation_observed2, fill_value, label_column,
                                plot_text_angle, plot_text_size, point_size,
                                gene_colour, sequence_data,
                                plot_sidechain=FALSE,layers=NULL,
                                paletteA=NULL, paletteB=NULL)
{
  
  # build the various features of the lolliplot
  
  # Build gene base either using domain information
  # or AA sidechain information
  if(plot_sidechain == TRUE)
  {
    sequence_data$coord_start <-
      as.numeric(as.character(sequence_data$coord)) - 1
    sequence_data$coord_end <- as.numeric(as.character(sequence_data$coord))
    gene_plot <- geom_rect(data=sequence_data,
                           mapping=aes_string(xmin='coord_start',
                                              xmax='coord_end', ymin=-.1,
                                              ymax=.1, fill='sidechain'))
    domain_plot <- geom_blank()
  } else {
    gene_plot <- geom_rect(data=gene_data[1,],
                           mapping=aes_string(xmin='pos_from',
                                              xmax='pos_to',
                                              ymin='height_min',
                                              ymax='height_max'),
                           fill='#999999', colour='#000000')
    # Take into account there might not be any domains
    if(nrow(gene_data) == 1)
    {
      domain_plot <- geom_blank()
    } else {
      domain_plot <- geom_rect(data=gene_data[-1,],
                               mapping=aes_string(xmin='pos_from',
                                                  xmax='pos_to',
                                                  ymin='height_min',
                                                  ymax='height_max',
                                                  fill='Domain'),
                               alpha=1, colour='black')
    }
  }
  
  # Build the Observed track
  observed_plot <- geom_point(data=mutation_observed,
                              mapping=aes_string(x='coord_x_dodge',
                                                 y='coord_y_dodge',
                                                 colour=fill_value),
                              size=point_size)
  observed_line <- geom_segment(data=mutation_observed,
                                mapping=aes_string(x='mutation_coord',
                                                   y=.1, xend='coord_x_dodge',
                                                   yend=.3))
  observed_line_2 <- geom_segment(data=mutation_observed,
                                  mapping=aes_string(x='coord_x_dodge', y=.3,
                                                     xend='coord_x_dodge',
                                                     yend='coord_y_dodge'))
  
  # Miscelaneous features
  title <- ggtitle(gene_data[1,1])
  x_label <- xlab('Amino Acid Position')
  
  # add a theme and guide to the plot
  theme <- theme(legend.position='bottom',
                 legend.direction='vertical',
                 legend.box='horizontal',
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title.y=element_blank())
  guide <- guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))
  
  # set colours manually if these are specified
  if(!is.null(paletteA))
  {
    gene_features_fill <- scale_fill_manual(values=paletteA)
  } else {
    gene_features_fill <- geom_blank()
  }
  
  if(!is.null(paletteB))
  {
    lolli_features_fill <- scale_colour_manual(values=paletteB)
  } else {
    lolli_features_fill <- geom_blank()
  }
  
  # construct the plot with or without 2nd observed track
  if(is.null(mutation_observed2))
  {
    y_limits <- ylim(c(-.1, max(mutation_observed$coord_y_dodge) + .1))
    tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
    p1 <- ggplot(data=tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
      gene_plot + domain_plot + observed_line_2 + observed_line +
      observed_plot + x_label + title + y_limits + theme_bw() + theme +
      guide + layers + gene_features_fill + lolli_features_fill
  } else {
    y_limits <- ylim(c(min(mutation_observed2$coord_y_dodge) - .1,
                       max(mutation_observed$coord_y_dodge) + .1))
    if(any(colnames(mutation_observed2) %in% fill_value))
    {
      observed2_plot <- geom_point(data=mutation_observed2,
                                   mapping=aes_string(x='coord_x_dodge',
                                                      y='coord_y_dodge',
                                                      colour=fill_value),
                                   size=point_size)
    } else {
      observed2_plot <- geom_point(data=mutation_observed2,
                                   mapping=aes_string(x='coord_x_dodge',
                                                      y='coord_y_dodge'),
                                   size=point_size)
    }
    observed2_line <- geom_segment(data=mutation_observed2,
                                   mapping=aes_string(x='mutation_coord',
                                                      y=-.1,
                                                      xend='coord_x_dodge',
                                                      yend=-.3))
    observed2_line_2 <- geom_segment(data=mutation_observed2,
                                     mapping=aes_string(x='coord_x_dodge',
                                                        y=-.3,
                                                        xend='coord_x_dodge',
                                                        yend='coord_y_dodge'))
    
    tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
    p1 <- ggplot(tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
      gene_plot + domain_plot + observed_line + observed_line_2 +
      observed_plot + observed2_line + observed2_line_2 + observed2_plot +
      x_label + title + y_limits + theme_bw() + theme + guide + layers +
      gene_features_fill + lolli_features_fill
  }
  
  # If a label column is specified plot labels
  if(any(colnames(mutation_observed) %in% "labels"))
  {
    mutation_observed$y_label_offset <-
      mutation_observed$coord_y_dodge + .01
    p1 <- p1 + geom_text(data=mutation_observed,
                         mapping=aes_string(x='coord_x_dodge',
                                            y='y_label_offset',
                                            label='labels'),
                         angle=plot_text_angle, size=plot_text_size,
                         vjust=1, hjust=0)
  }
  if(any(colnames(mutation_observed2) %in% "labels"))
  {
    mutation_observed2$y_label_offset <-
      mutation_observed2$coord_y_dodge - .01
    p1 <- p1 + geom_text(data=mutation_observed2,
                         mapping=aes_string(x='coord_x_dodge',
                                            y='y_label_offset',
                                            label='labels'),
                         angle=plot_text_angle, size=plot_text_size,
                         vjust=0, hjust=1)
  }
  
  return(p1)
}

#' Convert Codon to AA
#' 
#' Convert a Codon to the appropriate amino acid
#' @name lolliplot_Codon2AA
#' @param x Character string of length 1 giving the DNA codon to convert
#' @return Character corresponding to the residue for the given codon
#' @noRd

lolliplot_Codon2AA <- function(x)
{
  # Convert codons to single AA code
  x <- toupper(x)
  x <- switch(x, "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", "CTT"="L",
              "CTC"="L", "CTA"="L", "CTG"="L", "ATT"="I", "ATC"="I",
              "ATA"="I", "ATG"="M", "GTT"="V", "GTC"="V", "GTA"="V",
              "GTG"="V", "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
              "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P", "ACT"="T",
              "ACC"="T", "ACA"="T", "ACG"="T", "GCT"="A", "GCC"="A",
              "GCA"="A", "GCG"="A", "TAT"="Y", "TAC"="Y", "TAA"="OCHRE",
              "TAG"="AMBER", "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
              "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", "GAT"="D",
              "GAC"="D", "GAA"="E", "GAG"="E", "TGT"="C", "TGC"="C",
              "TGA"="OPAL", "TGG"="W", "CGT"="R", "CGC"="R", "CGA"="R",
              "CGG"="R", "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
              "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G")
  
  return(x)
}

#' Construct gene information
#' 
#' Build gene for input into lolliplot_buildMain
#' @name lolliplot_constructGene
#' @param gene character string specifying gene name
#' @param domain_data object of class data frame specifying protien domain
#' information, obtained from lolliplot_fetchDomain, should contain columns
#' giving "description", "start", "end"
#' @param length integer specifying length of transcript in amino acids
#' @return object of class data frame giving gene and domain information
#' @noRd

lolliplot_constructGene <- function(gene, domain_data, length)
{
  # message
  message("Constructing gene track")
  
  # Construct basic gene information, if there are no domains return the gene
  gene <- data.frame(Domain=gene, pos_from=1, pos_to=length, nest=1)
  if(nrow(na.omit(domain_data)) == 0)
  {
    gene$height_min <- .1/(as.numeric(gene$nest))
    gene$height_max <- -.1/(as.numeric(gene$nest))
    gene$pos_from <- as.numeric(gene$pos_from)
    gene$pos_to <- as.numeric(gene$pos_to)
    return(gene)
  }
  
  # rename columns for domain_data and make sure description column is
  # not a factor
  colnames(domain_data) <- c("description", "start", "end")
  domain_data$description <- as.character(domain_data$description)
  
  # quality check of domain data
  if(max(domain_data$end) > length)
  {
    memo <- paste0("The end position of a domain: ",  max(domain_data$end),
                   " is exceeding the length of the protein:", length)
    warning(memo)
  } else if(min(domain_data$start) < 1) {
    memo <- paste0("The start position of a domain:",
                   min(domain_data$start),
                   "is less than the start of the protein", 1)
    warning(memo)
  }
  
  # Check that start coordinates are always less than the end coordinates
  if(any(domain_data$start >= domain_data$end))
  {
    memo <- paste0("Found a start position greater than an end position",
                   " in the protein features track. Check input to Z or",
                   "results of the biomaRt query using dataOut==TRUE.")
    warning(memo)
  }
  
  # determine which regions are overlapping and annotate which nest domain is
  # sort on start
  domain_data$start <- as.numeric(domain_data$start)
  domain_data$end <- as.numeric(domain_data$end)
  domain_data <- domain_data[order(domain_data$start),]
  
  # annotate nests
  nest <- vector('numeric')
  end <- vector('numeric')
  for(i in 1:nrow(domain_data))
  {
    # Remove from end any values <= gene$start[i]
    idx <- domain_data$start[i] < end
    end <- end[idx]
    
    nest <- c(nest, length(end))
    end <- c(end, domain_data$end[i])
  }
  
  # add this nest information to the data frame
  domain_data$nest <- nest + 1
  colnames(domain_data) <- c("Domain", "pos_from", "pos_to", "nest")
  
  # combine gene and domain information
  gene <- rbind(gene, domain_data)
  
  # annotate display heights based on nesting and make sure coord are numeric
  gene$height_min <- .1/(as.numeric(gene$nest))
  gene$height_max <- -.1/(as.numeric(gene$nest))
  gene$pos_from <- as.numeric(gene$pos_from)
  gene$pos_to <- as.numeric(gene$pos_to)
  
  return(gene)
}

#' Convert DNA character string
#' 
#' Convert a character string of nucleotides to amino acids or side chain class
#' @name lolliplot_DNAconv
#' @param x Character string of nucleotides to convert
#' @param to Character string specifying conversion to do, one of "codon", 
#' "residue", "sidechain"
#' @return Converted string of nucleotides as character vector
#' @noRd

lolliplot_DNAconv <- function(x, to="residue")
{
  # check if given character string is a multiple of 3
  if(nchar(x)%%3 != 0)
  {
    memo <- paste0("Coding sequence retrieved for given ensembl transctipt",
                   ", is not a multiple of three. output may not be,",
                   " accurate!")
    warning(memo)
  }
  
  # split the character string into codons
  codon <- substring(x, seq(1,nchar(x), 3), seq(3, nchar(x), 3))
  if(toupper(to)=="CODON")
  {
    return(as.character(codon))
  }
  
  # convert the codons into amino acid residues
  residue <- sapply(codon, lolliplot_Codon2AA)
  if(toupper(to)=="RESIDUE")
  {
    return(as.character(residue))
  }
  
  # convert the residues into sidechain classifications
  sidechain <- sapply(residue, lolliplot_AA2sidechain)
  if(toupper(to)=="SIDECHAIN")
  {
    return(as.character(sidechain))
  }
  
  # return a warning if code gets this far
  memo <- paste0("did not recognize input to variable \"to\",",
                 " returning residue data")
  warning(memo)
  return(as.character(residue))
}

#' dodge coordinates
#'
#' given amino acid position dodge on x axis
#' @name lolliplot_dodgeCoordX
#' @param x numeric vector of position coordinates on x axis
#' @param rep.fact repulsive factor for plotted mutations observed track
#' @param rep.dist.lmt repulsive distance limit for plotted mutations observed
#' track
#' @param attr.fact attraction factor for plotted mutations observed track
#' @param adj.max maximum position change for each iteration observed track
#' @param adj.lmt position adjustment limit which simulation stops observed
#' track
#' @param iter.max maximum iterations beyond which to stop the simulation
#' observed track
#' @return numeric vector of dodged position coordinates on x axis
#' @importFrom FField FFieldPtRep
#' @noRd

lolliplot_dodgeCoordX <- function(x, rep.fact=5000, rep.dist.lmt=500,
                                  attr.fact=.1, adj.max=.1, adj.lmt=.5,
                                  iter.max=50000)
{
  # Format into data frame with columns as x and y
  x <- as.data.frame(cbind(x, 0))
  colnames(x) <- c('x', 'y')
  
  # Forcefield does not work with only 1 point, check for this first
  if(nrow(x) < 2)
  {
    return(x$x)
  }
  
  # take the data frame and apply a repulsive force to coordinates
  x <- FField::FFieldPtRep(x, rep.fact=rep.fact, rep.dist.lmt=rep.dist.lmt,
                           attr.fact=attr.fact, adj.max=adj.max,
                           adj.lmt=adj.lmt, iter.max=iter.max)
  
  return(x$x)
}

#' dodge coordinates
#'
#' given a data frame, dodge x coordinates ontop of each other
#' @name lolliplot_dodgeCoordY
#' @param x data frame containing columns coord_x_dodge
#' @param track character vector, one of "top", "bottom" specifying whether to
#' dodge in a positive or negative fashion
#' @return numeric vector of dodged position coordinates on y axis
#' @noRd

lolliplot_dodgeCoordY <- function(x, track='top')
{
  for(i in 1:length(x$coord_x_dodge))
  {
    if(track == 'top' & i == 1)
    {
      pos <- .3
      orig_pos <- .3
      pos_change <- .1
    } else if(track == 'bottom' & i == 1) {
      pos <- -.3
      orig_pos <- -.3
      pos_change <- -.1
    }
    
    if(i == 1)
    {
      y_axis_vec <- c(pos)
      next
    } else {
      x_coord_a <- x$coord_x_dodge[i-1]
      x_coord_b <- x$coord_x_dodge[i]
    }
    
    if(x_coord_b == x_coord_a)
    {
      new_pos <- pos + pos_change
      y_axis_pos <- new_pos
      y_axis_vec <- c(y_axis_vec, y_axis_pos)
      pos <- new_pos
      next
    } else {
      pos <- orig_pos
      y_axis_vec <- c(y_axis_vec, pos)
    }
  }
  
  return(y_axis_vec)
}

#' fetch protein domains
#' 
#' Retrieve protein domains given ensembl transcript ID
#' @name lolliplot_fetchDomain
#' @param transcriptID String specifying ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return data frame of protien domains and start/stop coordinates
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @noRd
lolliplot_fetchDomain <- function(transcriptID,
                                  species="hsapiens",
                                  host="www.ensembl.org")
{
  # display message
  message("Querying biomaRt for protein domains")
  
  # Load in mart
  ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                   host="www.ensembl.org")
  
  # select proper data set given regexp print warnings if unexpected out occur
  dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
  index <- which(grepl(species, dataset))
  if(length(index)>1)
  {
    memo <- paste0(species, " Matches more than one dataset for the",
                   " ensembl mart, please specify a species in the, ",
                   "following format: hsapiens")
    stop(memo)
  } else if(length(index)==0) {
    memo <- paste0(species, " does not appear to be supported by biomaRt",
                   "if you beleive this to be in error please modify", 
                   "you're input to to conform to this format: hsapiens")
    stop(memo)
  }
  ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                      mart=ensembl_mart)
  
  # Apply various filters using vector of values
  filters <- c("ensembl_transcript_id")
  values <- as.list(c(transcriptID))
  
  # Select attributes to retrieve (protein domain, start, stop)
  attributes <- c("interpro_description",
                  "interpro_start",
                  "interpro_end")
  
  # Retrieve data
  result <- biomaRt::getBM(attributes=attributes, filters=filters,
                           values=values, mart=ensembl_mart)
  
  return(result)
}

#' format mutation observations
#'
#' Create a data frame of mutation observations
#' @name lolliplot_mutationObs
#' @param x object of class data frame with columns trv_type and amino acid
#' change
#' @param track character string specifying one to 'top', 'bottom' to specify
#' proper track
#' @param fill_value character string giving the name of the column to shade
#' variants on
#' @param label_column character string specifying column containing text
#' information to be plotted
#' @param rep.fact repulsive factor for plotted mutations observed track
#' @param rep.dist.lmt repulsive distance limit for plotted mutations observed
#' track
#' @param attr.fact attraction factor for plotted mutations observed track
#' @param adj.max maximum position change for each iteration observed track
#' @param adj.lmt position adjustment limit which simulation stops observed
#' track
#' @param iter.max maximum iterations beyond which to stop the simulation
#' observed track
#' @return object of class data frame giving mutation observations
#' @noRd

lolliplot_mutationObs <- function(x, track, fill_value, label_column,
                                  rep.fact, rep.dist.lmt, attr.fact, adj.max,
                                  adj.lmt, iter.max)
{
  # Remove variants within an intronic or splice site region
  if(any(grepl("^e", x$amino_acid_change, ignore.case=TRUE, perl=TRUE)))
  {   
    # save original data frame size before subset for message
    origDim <- nrow(x)
    
    # remove regions with AA change starting with e (i.e. intronic/splice)
    x <- x[-which(grepl("^e", x$amino_acid_change,
                        ignore.case=TRUE, perl=TRUE)),]
    
    newDim <- nrow(x)
    
    # print update message
    memo <- paste0("Removed ", origDim - newDim,
                   " variants not within a residue")
    message(memo)
    
    # if the removal has removed all rows print an error
    if(newDim == 0)
    {
      memo <- paste0("Did not detect any residues, please check input", 
                     " lolliplot must have at least one valid residue",
                     " present!")
      stop(memo)
    }
  }
  
  # extract the mutation types and set a flag specifying they are present
  if(any(colnames(x) %in% fill_value))
  {
    fill_value_flag <- TRUE
    fill <- as.character(x[,eval(fill_value)])
  } else {
    fill_value_flag <- FALSE
  }
  
  # extract the mutation coordinates
  mutation_coord <- x$amino_acid_change
  if(all(grepl("p\\.", mutation_coord)))
  {
    message("Detected p. notation for amino_acid_change")
    mutation_coord <- as.numeric(gsub("p\\.[*a-zA-z]*(\\d+).*?$", "\\1",
                                      mutation_coord, perl=TRUE))
  } else if(all(grepl("c\\.", mutation_coord))) {
    memo <- paste0("c. notation is not currently supported",
                   " please specify amino acid change in p. notation")
    stop(memo)
  } else {
    memo <- paste0("Could not determine notation type for ",
                   "column \"amino_acid_change\", please check input.", 
                   "Expecting p. notation: ex. p.R383A")
    stop(memo)
  }
  
  # combine mutation type and mutation coord into a data frame
  if(fill_value_flag)
  {
    mutation_data <- as.data.frame(cbind(mutation_coord, fill))
    colnames(mutation_data) <- c('mutation_coord', eval(fill_value))
  } else {
    mutation_data <- as.data.frame(mutation_coord)
    colnames(mutation_data) <- c('mutation_coord')
  }
  mutation_data$mutation_coord <-
    as.numeric(as.character(mutation_data$mutation_coord))
  
  # add extra column giving height of Y axis for points to be plotted
  if(track == 'top')
  {
    mutation_data$height_max <- .3
  } else if (track == 'bottom') {
    mutation_data$height_min <- -.3
  } else {
    stop("Fatal error: incorrect track type specified")
  }
  
  # extract the mutation types and set a flag specifying they are present
  if(any(colnames(x) %in% label_column))
  {
    label_column_flag <- TRUE
    mutation_data$labels <- as.character(x[,eval(label_column)])
  } else {
    label_column_flag <- FALSE
  }
  
  # Dodge mutation coordinates on the x axis
  if(track == 'top')
  {
    memo <- paste0("applying force field to observed mutations for",
                   " top track. This will take time if n is large",
                   ", see vignette for tips")
    message(memo)
  } else if (track == 'bottom') {
    memo <- paste0("applying force field to observed mutations for",
                   " bottom track. This will take time if n is large",
                   ", see vignette for tips")        
    message(memo)
  }
  mutation_data <- mutation_data[order(mutation_coord),] 
  mutation_data$coord_x_dodge <- 
    lolliplot_dodgeCoordX(as.vector(mutation_data$mutation_coord),
                          rep.fact=rep.fact,
                          rep.dist.lmt=rep.dist.lmt,
                          attr.fact=attr.fact,
                          adj.max=adj.max,
                          adj.lmt=adj.lmt,
                          iter.max=iter.max)
  
  # Dodge y coordinates
  if(track == 'top')
  {
    mutation_data$coord_y_dodge <- lolliplot_dodgeCoordY(mutation_data,
                                                         track='top')
  } else if(track == 'bottom') {
    mutation_data$coord_y_dodge <- lolliplot_dodgeCoordY(mutation_data,
                                                         track='bottom')
  }
  
  return(mutation_data)
}

#' Check input to lolliplot
#'
#' Perform Basic quality checks for lolliplot input
#' @name lolliplot_qual
#' @param x object of class data frame containing columns transcript_name, gene,
#' and amino_acid_change and rows denoting mutations
#' @param y object of class data frame containing columns transcript_name, and
#' amino_acid_change and rows denoting mutations
#' @param z Object of class data frame containing columns "description", "start",
#'  "stop" specifying gene regions to highlight
#' @return objects passing basic quality checks
#' @noRd

lolliplot_qual <- function(x, y, z)
{
  # Check input to x
  if(!is.data.frame(x))
  {
    message("Input to x is not a data frame, attempting to coerce")
    x <- as.data.frame(x)
    x <- droplevels(x)
  }
  
  # Check for correct columns in x
  if(!all(c('transcript_name', 'gene', 'amino_acid_change') %in% colnames(x)))
  {
    stop("Did not detect correct columns in x,
         missing one of transcript_name, gene, amino_acid_change")
  }
  
  # Make sure columns in x are of the proper class
  x$transcript_name <- as.factor(x$transcript_name)
  x$gene <- as.factor(x$gene)
  x$amino_acid_change <- as.factor(x$amino_acid_change)
  
  # Check that "transcript_name" in x contains only 1 transcript
  if(length(unique(x$transcript_name)) != 1)
  {
    stop("Detected more than 1 transcript in input to x")
  }
  
  # Check input to y
  if(!is.null(y))
  {
    # is y a data frame?
    if(!is.data.frame(y))
    {
      message(y, "is not a data frame, attempting to coerce")
      y <- as.data.frame(y)
      y <- droplevels(y)
    }
    
    # does y have correct columns?
    if(!all(c('transcript_name', 'amino_acid_change') %in% colnames(y)))
    {
      stop("Did not detect correct columns in y, missing one of
           transcript_name, amino_acid_change")
    }
    
    # make sure columns in y are of proper class
    y$transcript_name <- as.factor(y$transcript_name)
    y$amino_acid_change <- as.factor(y$amino_acid_change)
    }
  
  # Check input to z
  if(!is.null(z))
  {
    # is z a data frame?
    if(!is.data.frame(z))
    {
      memo <- paste0("Input to z is not a data frame",
                     ", attempting to coerce")
      warning(memo)
      z <- as.data.frame(z)
      z <- droplevels(z)
    }
    
    # does z contain correct columns
    if(!all(c("description", "start", "stop") %in% colnames(z)))
    {
      memo <- paste0("Did not detect correct columns in input to z, ",
                     "missing one of \"description\", \"start\",",
                     " \"stop\"")
      stop(memo)
    }
    
    # make sure column class of z are of proper type
    z$description <- as.factor(z$description)
    z$start <- as.numeric(as.character(z$start))
    z$stop <- as.numeric(as.character(z$stop))
  }
  
  return(list(x, y, z))
  }

#' Reduce Lolli
#' 
#' Reduce lollis stacked ontop of each other to the amount specified
#' @name lolliplot_reduceLolli
#' @param x Data frame with column name mutation_coord to reduce lollis on
#' @param max Integer specifying the maximum number of lollis to allow
#' @return Object of class data frame taking the reduced form of x
#' @importFrom plyr count
#' @noRd

lolliplot_reduceLolli <- function(x, max=NULL)
{
  # if max is null no reduction is required
  if(is.null(max))
  {
    return(x)
  }
  
  # Get a frequency of counts from x
  coordFreq <- plyr::count(x$mutation_coord)
  colnames(coordFreq) <- c("x", "freq")
  
  keep <- vector('numeric')
  # Loop through mutations keeping only what is under max
  for(i in 1:nrow(coordFreq))
  {
    index <- which(x$mutation_coord == coordFreq$x[i])
    index <- index[1:max]
    keep <- c(keep, index)
  }
  
  # subset the input on what we want to keep
  x <- x[keep,]
  
  return(x)
}

#' fetch protein length
#' 
#' Retrieve protein length from ensembl database given enseml transcript id
#' @name lolliplot_transcriptID2codingSeq
#' @param transcriptID character string giving ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return length in residues of ensembl transcript id
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @noRd

lolliplot_transcriptID2codingSeq <- function(transcriptID,
                                             species="hsapiens",
                                             host="www.ensembl.org")
{
  # display mesage
  memo <- paste0("Using the following host: ", host, " for biomaRt queries",
                 " to change the ensembl annotation version alter this",
                 " parameter!")
  message(memo)
  message("Querying biomaRt for transcript sequence")
  
  # Load in mart
  ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                   host=host)
  
  # select proper data set given regexp print warnings if unexpected out occur
  dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
  index <- which(grepl(species, dataset))
  if(length(index)>1)
  {
    memo <- paste0(species, " Matches more than one dataset for the",
                   " ensembl mart, please specify a species in the, ",
                   "following format: hsapiens")
    stop(memo)
  } else if(length(index)==0) {
    valid_species <- toString(gsub("_gene_ensembl",
                                   "",
                                   dataset))
    
    memo <- paste0(species, " does not appear to be supported by biomaRt",
                   " please specify one of the following species:",
                   valid_species)
    stop(memo)
  }
  ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                      mart=ensembl_mart)
  
  # Apply various filters using vector of values
  filters <- c("ensembl_transcript_id")
  ensg_id <- as.character(transcriptID)
  
  # Select attributes to retrieve coding dna sequence
  attributes <- as.list(c("coding","cds_length"))
  
  # Retrieve data
  result <- biomaRt::getBM(attributes=attributes, filters=filters,
                           values=ensg_id, mart=ensembl_mart)
  
  return(as.list(result))
}

#' Choose output
#' 
#' Selector for choosing output for GenVisR functions
#' 
#' @name multi_selectOut
#' @param data Data object to output
#' @param plot Plot object to output
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot".
#' @param draw Boolean specifying if the input to plot needs to be drawn
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid grid.draw
#' @noRd

multi_selectOut <- function(data, plot, out="plot", draw="FALSE")
{
  # Decide what to output
  if(toupper(out) == "DATA")
  {
    return(data)
  } else if(toupper(out) == "PLOT" & isTRUE(draw)) {
    return(grid::grid.draw(plot))
  } else if(toupper(out) == "PLOT" & !isTRUE(draw)) {
    return(plot)
  } else if(toupper(out) == "GROB" & isTRUE(draw)) {
    return(plot)
  } else if(toupper(out) == "GROB" & !isTRUE(draw)) {
    return(ggplot2::ggplotGrob(plot))
  } else {
    warning("Did not recognize input to out...")
    if(isTRUE(draw))
    {
      return(grid::grid.draw(plot))
    } else {
      return(plot)
    }
  }
}

