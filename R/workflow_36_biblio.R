# ==========================================================================
# workflow of biblio
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_biblio <- setClass("job_biblio", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/massimoaria/bibliometrix"),
    cite = "[@BibliometrixAria2017]",
    method = "Package bibliometrix used for quantitative research in bibliometrics and scientometrics"
    ))

job_biblio <- function(lst)
{
  kds <- do.call(wos_bind, lapply(kds, wos_or))
  x <- .job_biblio()
  x$keywords <- kds
  return(x)
}

setMethod("step0", signature = c(x = "job_biblio"),
  function(x){
    step_message("Prepare your data with function `job_biblio`.")
  })

setMethod("step1", signature = c(x = "job_biblio"),
  function(x, num, from = "~/Downloads", to = timeName("wos"), pattern = "^savedrecs"){
    step_message("Re-save the obtained files from WOS.")
    files <- unlist(moveToDir(paste0("wos_", 1:num), pattern, from = "~/Downloads",
      to = to, suffix = ".txt"))
    object(x) <- e(bibliometrix::convert2df(files, dbsource = "wos", format = "plaintext"))
    alls <- as_tibble(data.frame(object(x)))
    x$summary <- e(bibliometrix::missingData(object(x)))
    x@tables[[ 1 ]] <- namel(alls)
    x$files <- files
    x$savedir <- to
    return(x)
  })

setMethod("step2", signature = c(x = "job_biblio"),
  function(x){
    step_message("Perform a descriptive analysis of the bibliographic data frame.")
    require(bibliometrix)
    x$overview <- e(bibliometrix::biblioAnalysis(object(x)))
    cli::cli_alert_info("bibliometrix:::plot.bibliometrix")
    p.most_ave <- .plot_biblioAnalysis(x$overview)
    x@plots[[ 2 ]] <- namel(p.most_ave)
    return(x)
  })

setMethod("step3", signature = c(x = "job_biblio"),
  function(x, n = 30){
    step_message("Network analysis.")
    require(bibliometrix)
    ## country
    if (T) {
      obj <- e(bibliometrix::metaTagExtraction(object(x), Field = "AU_CO"))
      obj.net <- e(bibliometrix::biblioNetwork(obj, analysis = "collaboration",
          network = "countries"))
      p.net <- e(bibliometrix::networkPlot(obj.net, n = dim(obj.net)[1],
          Title = "Country Collaboration", type = "circle", size = TRUE,
          remove.multiple = FALSE, labelsize = 0.8))
      p.co_country <- .plot_coNetwork(p.net$graph, fill = "color", size = "deg") +
        guides(fill = "none", edge_width = "none") +
        labs(size = "Centrality degree")
      p.co_country <- wrap(p.co_country, 9, 7)
    }
    ## citation
    if (T) {
      ## C/Y filter and export
      obj.net <- e(bibliometrix::biblioNetwork(object(x), analysis = "co-citation",
          network = "references", n = n))
      ## the most cited
      rownames(obj.net)
      strsplit(x@tables$step1$alls$rownames, ", ")

      p.net <- e(bibliometrix::networkPlot(obj.net,
          Title = "Co-Citation Network", type = "fruchterman",
          size = T, remove.multiple = FALSE, labelsize = 0.7, edgesize = 5))
      p.co_citation <- .plot_coNetwork(p.net$graph, fill = "color", size = "deg",
        edge_color = "orange", zoRange = 1.5, fun_edge = ggraph::geom_edge_fan,
        width_range = c(.01, .05)
      )
      p.co_citation <- p.co_citation +
        guides(fill = "none", edge_width = "none") +
        labs(size = "Centrality degree")
      p.co_citation <- wrap(p.co_citation, 9, 7)
    }
    ## keywords
    obj.net <- e(bibliometrix::biblioNetwork(object(x), analysis = "co-occurrences", network = "keywords"))
    p.net <- e(bibliometrix::networkPlot(obj.net, normalize = "association", weighted = T, n = n,
      Title = "Keyword Co-occurrences", type = "fruchterman", size = T,
      edgesize = 5, labelsize = 0.7))
    p.co_keywords <- recordPlot()
    x@plots[[ 3 ]] <- namel(p.co_country, p.co_citation, p.co_keywords)
    return(x)
  })

setMethod("step4", signature = c(x = "job_biblio"),
  function(x){
    step_message("Co-Word Analysis and Historical Direct Citation Network.")
    # ID: Keywords Plus associated by ISI or SCOPUS database
    require(bibliometrix)
    cs <- e(bibliometrix::conceptualStructure(object(x), field = "ID", method = "MCA",
        minDegree = 10, clust = 5, stemming = FALSE, labelsize = 3,
        documents = 20, graph = FALSE))
    p.conceptual_structure <- cs$graph_terms +
      theme(
        axis.title = element_text(size = 15, face = "plain"),
        plot.title = element_text(size = 15, face = "plain")
      )
    p.conceptual_structure$layers[[ 6 ]] <- NULL
    p.conceptual_structure <- wrap(p.conceptual_structure, 12, 10)
    histResults <- e(bibliometrix::histNetwork(object(x)))
    p.hist <- e(bibliometrix::histPlot(histResults, n = 10, size = FALSE, label = "short"))
    p.hist$g$layers[[ 4 ]] <- NULL
    p.hist$g <- p.hist$g + coord_cartesian(xlim = zoRange(p.hist$g$data$x, 1.2))
    p.hist_direct <- wrap(p.hist$g)
    x$cs <- cs
    x$hist_direct <- p.hist
    x@plots[[ 4 ]] <- namel(p.conceptual_structure, p.hist_direct)
    return(x)
  })

.plot_coNetwork <- function(igraph, label = "name",
  fill = "name", size = "size",
  width = "width", edge_color = "lightblue", zoRange = 1.2,
  fun_edge = ggraph::geom_edge_arc, width_range = c(.5, 1))
{
  graph <- tidygraph::as_tbl_graph(igraph)
  graph <- create_layout(graph, layout = 'linear', circular = T)
  p <- ggraph(graph) +
    geom_node_point(aes(size = !!rlang::sym(size), fill = !!rlang::sym(fill)), alpha = .7, shape = 21) +
    geom_node_text(aes(x = x * 1.1,  y = y * 1.1, 
        label = !!rlang::sym(label),
        angle = -((-node_angle(x,  y) + 90) %% 180) + 90), 
      size = 3, hjust = 'outward') +
    fun_edge(aes(width = !!rlang::sym(width)), color = edge_color, alpha = .5) +
    scale_size(range = c(3, 9)) +
    scale_fill_manual(values = color_set()) +
    scale_edge_width_continuous(range = width_range) +
    coord_cartesian(xlim = zoRange(graph$x, zoRange), ylim = zoRange(graph$y, zoRange)) +
    theme_graph()
  p
}

wos_logic <- function(lst, f = "ALL", logic = c("AND", "OR", "NOT")) {
  logic <- match.arg(logic)
  kds <- unlist(lst)
  kds <- paste0(f, "=(", kds, ")")
  res <- kds[1]
  for (i in 2:length(kds)) {
    res <- paste0("(", res, ")", " ", logic, " ", kds[i])
  }
  res
}

wos_and <- function(lst, f = "ALL") {
  wos_logic(lst, f, logic = "AND")
}

wos_or <- function(lst, f = "ALL") {
  wos_logic(lst, f, logic = "OR")
}

wos_not <- function(lst, f = "ALL") {
  wos_logic(lst, f, logic = "NOT")
}

wos_bind <- function(..., logic = c("AND", "OR", "NOT")) {
  logic <- match.arg(logic)
  ct <- unlist(list(...))
  res <- ct[1]
  for (i in 2:length(ct)) {
    res <- paste0("(", res, ")", " ", logic, " (", ct[i], ")")
  }
  res
}


## the following copied from other packages
.plot_biblioAnalysis <- function(x, ...){
  graphs <- list()
  arguments <- list(...)
  if (sum(names(arguments)=="k")==0){k=10} else {k=arguments$k}
  if (sum(names(arguments)=="pause")==0){pause=FALSE} else {pause=arguments$pause}
  
  if (pause == TRUE){
    cat("Hit <Return> to see next plot: ")
    line <- readline()}
  
  xx=as.data.frame(x$Authors[1:k])
  xcoord <- c(k-0.2-(k)*0.15, k-0.02)+1
  ycoord <- c(max(xx$Freq),max(xx$Freq)-base::diff(range(xx$Freq))*0.15)
  # Authors
  
  g <- ggplot(data = xx, aes(x = reorder(AU, Freq), y = Freq)) +
    geom_bar(stat="identity", fill="grey90")+
    labs(title="Most productive Authors", x = "Authors")+
    labs(y = "N. of Documents")+
    theme(text = element_text(color = "#444444")
          ,panel.background = element_rect(fill = '#FFFFFF')
          ,panel.grid.minor = element_line(color = '#EFEFEF')
          ,panel.grid.major = element_line(color = '#EFEFEF')
          ,plot.title = element_text(size = 18)
          ,axis.title = element_text(size = 14, color = '#555555')
          ,axis.title.y = element_text(vjust = 1, angle = 90)
          ,axis.title.x = element_text(hjust = 0)
          ,axis.text.x = element_text(size=10)
          ,axis.line.x = element_line(color="black",size=0.5)
          ,axis.line.y = element_line(color="black",size=0.5)
    ) +
    coord_flip()
  
  graphs$MostProdAuthors <- g
  
  if (pause == TRUE){
    cat("Hit <Return> to see next plot: ")
    line <- readline()}
  if (!is.na(x$CountryCollaboration[1,1])){
  # Countries
  xx=x$CountryCollaboration[1:k,]
  xx=xx[order(-(xx$SCP+xx$MCP)),]
  xx1=cbind(xx[,1:2],rep("SCP",k))
  names(xx1)=c("Country","Freq","Collaboration")
  xx2=cbind(xx[,c(1,3)],rep("MCP",k))
  names(xx2)=c("Country","Freq","Collaboration")
  xx=rbind(xx2,xx1)
  xx$Country=factor(xx$Country,levels=xx$Country[1:dim(xx2)[1]])
  
  Freq <- x$CountryCollaboration$SCP[1:k]+x$CountryCollaboration$MCP[1:k]
  st <- floor(k/10)
  #xcoord <- c(st-0.2-(st)*0.85, 0.02)+1
  xcoord <- c(1,max(st,3))
  ycoord <- c(max(Freq),max(Freq)-base::diff(range(Freq))*0.15)

  g <- suppressWarnings(ggplot(data=xx, aes(x=.data$Country, y=.data$Freq,fill=.data$Collaboration)) +
    geom_bar(stat="identity")+
    scale_x_discrete(limits = rev(levels(xx$Country)))+
    scale_fill_discrete(name="Collaboration",
                        breaks=c("SCP","MCP"))+
    labs(title = "Most Productive Countries", x = "Countries", y = "N. of Documents", 
         caption = "SCP: Single Country Publications, MCP: Multiple Country Publications")+
      theme(text = element_text(color = "#444444")
            ,panel.background = element_rect(fill = '#FFFFFF')
            ,panel.grid.minor = element_line(color = '#EFEFEF')
            ,panel.grid.major = element_line(color = '#EFEFEF')
            ,plot.title = element_text(size = 18)
            ,axis.title = element_text(size = 14, color = '#555555')
            ,axis.title.y = element_text(vjust = 1, angle = 90)
            ,axis.title.x = element_text(hjust = 0)
            ,axis.text.x = element_text(size=10)
            ,axis.line.x = element_line(color="black",size=0.5)
            ,axis.line.y = element_line(color="black",size=0.5)
      ) +
      coord_flip()) 
  graphs$MostProdCountries <- g
  } else {
    graphs$MostProdCountries <- NA
  }
  
  if (pause == TRUE){
    cat("Hit <Return> to see next plot: ")
    line <- readline()
  }

  # Articles per Year
  
  Tab=table(x$Years)
  
  ## inserting missing years
  YY=setdiff(seq(min(x$Years, na.rm=TRUE),max(x$Years, na.rm=TRUE)),names(Tab))
  Y=data.frame(Year=as.numeric(c(names(Tab),YY)),Freq=c(as.numeric(Tab),rep(0,length(YY))))
  Y=Y[order(Y$Year),]
  
  names(Y)=c("Year","Freq")
  
  xcoord <- c(max(Y$Year)-0.02-base::diff(range(Y$Year))*0.15, max(Y$Year)-0.02)+1
  ycoord <- c(min(Y$Freq),min(Y$Freq)+base::diff(range(Y$Freq))*0.15)
  
  g <- ggplot(Y, aes(x = .data$Year, y = .data$Freq)) +
    geom_line() +
    geom_area(fill = 'grey90', alpha = .5) +
    labs(x = 'Year'
         , y = 'Articles'
         , title = "Annual Scientific Production") +
    scale_x_continuous(breaks= (Y$Year[seq(1,length(Y$Year),by=2)])) +
    theme(text = element_text(color = "#444444")
          ,panel.background = element_rect(fill = '#FFFFFF')
          ,panel.grid.minor = element_line(color = '#EFEFEF')
          ,panel.grid.major = element_line(color = '#EFEFEF')
          ,plot.title = element_text(size = 18)
          ,axis.title = element_text(size = 14, color = '#555555')
          ,axis.title.y = element_text(vjust = 1, angle = 90)
          ,axis.title.x = element_text(hjust = 0)
          ,axis.text.x = element_text(size=10, angle = 90)
          ,axis.line.x = element_line(color="black",size=0.5)
          ,axis.line.y = element_line(color="black",size=0.5)
    )
  
  graphs$AnnualScientProd <- g
  
  Table2=NA
  if(!(x$DB %in% c("COCHRANE","PUBMED"))){
  
  if (pause == TRUE){
    cat("Hit <Return> to see next plot: ")
    line <- readline()}
  
  # Total Citation Plot
  Table2=aggregate(x$TotalCitation,by=list(x$Years),length)
  Table2$xx=aggregate(x$TotalCitation,by=list(x$Years),mean)$x
  Table2$Annual=NA
  d=date()
  d=as.numeric(substring(d,nchar(d)-3,nchar(d)))
  Table2$Years=d-Table2$Group.1
  Table2$Annual=Table2$xx/Table2$Years
  names(Table2)=c("Year","N","MeanTCperArt","MeanTCperYear","CitableYears")
  
  ## inserting missing years
  YY=setdiff(seq(min(x$Years,na.rm=TRUE),max(x$Years,na.rm=TRUE)),Table2$Year)
  if (length(YY>0)){
    YY=data.frame(YY,0,0,0,0)
    names(YY)=c("Year","N","MeanTCperArt","MeanTCperYear","CitableYears")
    Table2=rbind(Table2,YY)
    Table2=Table2[order(Table2$Year),]
    row.names(Table2)=Table2$Year
  }

  xcoord <- c(max(Table2$Year)-0.02-base::diff(range(Table2$Year))*0.15, max(Table2$Year)-0.02)+1
  ycoord <- c(min(Table2$MeanTCperYear),min(Table2$MeanTCperYear)+base::diff(range(Table2$MeanTCperYear))*0.15)

  g <- ggplot(Table2, aes(x = .data$Year, y = .data$MeanTCperYear)) +
    geom_line() +
    geom_area(fill = 'grey90', alpha = .5) +
    labs(x = 'Year'
      , y = 'Citations'
      , title = "Average Article Citations per Year")+
    scale_x_continuous(breaks= (Table2$Year[seq(1,length(Table2$Year),by=2)])) +
    theme(text = element_text(color = "#444444")
      ,panel.background = element_rect(fill = '#FFFFFF')
      ,panel.grid.minor = element_line(color = '#EFEFEF')
      ,panel.grid.major = element_line(color = '#EFEFEF')
      ,plot.title = element_text(size = 18)
      ,axis.title = element_text(size = 14, color = '#555555')
      ,axis.title.y = element_text(vjust = 1, angle = 90)
      ,axis.title.x = element_text(hjust = 0)
      ,axis.text.x = element_text(size=10, angle = 90)
      ,axis.line.x = element_line(color="black",size=0.5)
      ,axis.line.y = element_line(color="black",size=0.5)
    )

    graphs$AverArtCitperYear <- g

    if (pause == TRUE){
      cat("Hit <Return> to see next plot: ")
      line <- readline()}

    xcoord <- c(max(Table2$Year)-0.02-base::diff(range(Table2$Year))*0.15, max(Table2$Year)-0.02)+1
    ycoord <- c(min(Table2$MeanTCperArt),min(Table2$MeanTCperArt)+base::diff(range(Table2$MeanTCperArt))*0.15)

    g <- ggplot(Table2, aes(x = .data$Year, y = .data$MeanTCperArt)) +
      geom_line() +
      geom_area(fill = 'grey90', alpha = .5) +
      labs(x = 'Year'
        , y = 'Citations'
        , title = "Average Total Citations per Year")+
      scale_x_continuous(breaks= (Table2$Year[seq(1,length(Table2$Year),by=2)])) +
      theme(text = element_text(color = "#444444")
        ,panel.background = element_rect(fill = '#FFFFFF')
        ,panel.grid.minor = element_line(color = '#EFEFEF')
        ,panel.grid.major = element_line(color = '#EFEFEF')
        ,plot.title = element_text(size = 18)
        ,axis.title = element_text(size = 14, color = '#555555')
        ,axis.title.y = element_text(vjust = 1, angle = 90)
        ,axis.title.x = element_text(hjust = 0)
        ,axis.text.x = element_text(size=10, angle = 90)
        ,axis.line.x = element_line(color="black",size=0.5)
        ,axis.line.y = element_line(color="black",size=0.5)
      )
      graphs$AverTotCitperYear <- g
  } else {
    graphs$AverArtCitperYear=NA
    graphs$AverTotCitperYear=NA
  }
  graphs
}

