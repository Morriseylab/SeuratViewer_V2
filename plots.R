require(dplyr)
monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


getGeneExpression <- function(cds,genes){
  markers_fData <- subset(fData(cds), gene_short_name %in% 
                            genes)
  cds_subset <- cds[row.names(markers_fData), ]
  #return(t(as.matrix(exprs(cds_subset))))  
  
  ### Adj values #####
  cds_exprs <- exprs(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
  cds_exprs <- as.matrix(log2(t(cds_exprs)))
  #cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  return(cds_exprs)
  
}


visualize_gene_markers2 <- function (cds, gene_probes, x=1,y=2, limits = c(0, 10), marker_size = 0.1, color='red',
                                     title = NULL) 
{
  
  gene_values<-getGeneExpression(cds,gene_probes)
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  
  projection <- data.frame(t(cds@reducedDimA[c(x, y), ]))
  
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- reshape2::melt(proj_gene, id.vars = c("Component.1", 
                                                          "Component.2"))
  
  
  
  p <- ggplot(arrange(proj_gene_melt, value), aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    facet_wrap(~variable) + scale_colour_gradient(low = "grey", 
                                                  high = color, name = "")
  
  
  p <- p + monocle_theme_opts() + xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) + theme(legend.position = "top", 
                                        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(text = element_text(size = 15))
  
  
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  # p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
  #                             panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

#################################
# Function to return max and min of gene to used to make silder
#################################
getGeneRange <- function(cds,gene_probes){
  gene_values<-getGeneExpression(cds,gene_probes)
  minr<- round(min(gene_values),2) 
  maxr<- round(max(gene_values),2)
  
  return(c(ifelse(minr==0,.1,minr-.1),maxr))
  #return(c(1,12))

}




bigene_getValues <- function(cds,gene_probes,limita,limitb){
  gene_values<-getGeneExpression(cds,gene_probes)
  colnames(gene_values) <- c('genea','geneb')
  
  as.data.frame(gene_values) %>% 
    mutate(value = ifelse(genea>=limita[1] & geneb>=limitb[1],
                          'both',
                          ifelse(genea>=limita[1] & geneb<limitb[1],
                                 gene_probes[1],
                                 ifelse(genea<=limita[1] & geneb>=limitb[1],
                                        gene_probes[2],
                                        'none')))
    ) #%>% select(value)

}




bigene_plot <- function (cds, gene_probes, x=1,y=2, limita=c(1,100), limitb=c(1,100), marker_size = 0.1, 
          title = NULL) 
{

  gene_values <- bigene_getValues(cds,gene_probes,limita,limitb)
  
  projection <- data.frame(t(cds@reducedDimA[c(x, y), ]))
  
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))

 proj_gene$value = factor(proj_gene$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
  
 proj_gene <- arrange(proj_gene, desc(value))
 
  p <- ggplot(proj_gene, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A", 'grey90'),drop=F) +
    theme(legend.key.size = unit(10,"point")) + xlab(paste("Component", x)) + 
    ylab(paste("Component", y))
 
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + monocle_theme_opts() + theme(plot.title = element_text(hjust = 0.5), 
                              legend.position="bottom", 
                              legend.title=element_blank(),
                              legend.text=element_text(size=14),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank())
  return(p)
}



bigene_pieChart <- function (gbm, gene_probes, limita, limitb, title = NULL)
{
  
  gene_values <- bigene_getValues(gbm,gene_probes,limita,limitb) %>% group_by(value) %>% dplyr::summarise(count=n()) 
  gene_values$value <-factor(gene_values$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
  gene_values<-arrange(gene_values, desc(value))
  pct <- round(gene_values$count/sum(gene_values$count)*100)
  lbls <- paste(gene_values$value,' (',pct, '%)',sep='') # add percents to labels 
  with(gene_values,pie(count,labels = lbls, col=rev(c("#E41A1C","#377EB8","#4DAF4A", 'grey90')),
      main=""))
  
}


