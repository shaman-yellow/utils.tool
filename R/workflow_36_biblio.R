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
    method = "R package `bibliometrix` used for quantitative research in bibliometrics and scientometrics",
    tag = "lite:bib",
    analysis = "Bibliometrix 文献计量分析"
    ))

job_biblio <- function(lst)
{
  kds <- do.call(wos_bind, lapply(lst, wos_or))
  x <- .job_biblio()
  keywords <- new_lich(c(lst, list(query = kds)))
  keywords <- .set_lab(keywords, "keywords", "for query")
  x$keywords <- keywords
  x$kds <- lst
  return(x)
}

setMethod("step0", signature = c(x = "job_biblio"),
  function(x){
    step_message("Prepare your data with function `job_biblio`.")
  })

setMethod("step1", signature = c(x = "job_biblio"),
  function(x, num, from = "~/Downloads", to = timeName("wos"), pattern = "^savedrecs",
    filter_dup = TRUE, filter_nonEng = TRUE, filter_naAbs = TRUE, filter_AbsKey = TRUE)
  {
    step_message("Re-save the obtained files from WOS.")
    files <- unlist(collateFiles(paste0("wos_", 1:num), pattern, from = "~/Downloads",
      to = to, suffix = ".txt"))
    object(x) <- e(bibliometrix::convert2df(files, dbsource = "wos", format = "plaintext"))
    dat <- alls <- as_tibble(data.frame(object(x)))
    dat <- dplyr::mutate(dat, .seq_id = seq_len(nrow(dat)))
    fun <- function(dat) {
      dplyr::filter(dat, !.seq_id %in% !!d.out$.seq_id)
    }
    if (filter_dup) {
      d.out <- dplyr::filter(dat, duplicated(TI))
      filter_dup <- nrow(d.out)
      dat <- fun(dat)
    }
    if (filter_nonEng) {
      d.out <- dplyr::filter(dat, LA != "ENGLISH")
      filter_nonEng <- nrow(d.out)
      dat <- fun(dat)
    }
    if (filter_naAbs) {
      d.out <- dplyr::filter(dat, is.na(AB))
      filter_naAbs <- nrow(d.out)
      dat <- fun(dat)
    }
    if (filter_AbsKey) {
      kds <- x$kds
      isThat.lst <- lapply(kds,
        function(ch) {
          grpl(dat$AB, paste0(ch, collapse = "|"), ignore.case = TRUE)
        })
      isThat <- isThat.lst[[1]]
      for (i in 2:length(isThat.lst)) {
        isThat <- isThat & isThat.lst[[ i ]]
      }
      notThat <- !isThat
      d.out <- dplyr::filter(dat, !!notThat)
      filter_AbsKey <- nrow(d.out)
      dat <- fun(dat)
    }
    keep <- dat$.seq_id
    object(x) <- object(x)[keep, ]
    x$summary <- e(bibliometrix::missingData(object(x)))
    x@tables[[ 1 ]] <- namel(alls)
    x$files <- files
    x$savedir <- to
    x$prefilter <- namel(filter_dup, filter_nonEng, filter_naAbs, filter_AbsKey)
    x$prefilter.lich <- new_lich(list(
      `Step 1: filter out the duplicated by Title` = filter_dup,
      `Step 2: filter out the non-Englich documents` = filter_nonEng,
      `Step 3: filter out the abstract of which documents is not available` = filter_naAbs,
      `Step 4: filter out the abstract of which not contain the keywords` = filter_AbsKey,
      `Finaly documents number` = nrow(object(x))
    ))
    return(x)
  })

setMethod("step2", signature = c(x = "job_biblio"),
  function(x, n = 10){
    step_message("Perform a descriptive analysis of the bibliographic data frame.")
    require(bibliometrix)
    x$overview <- e(bibliometrix::biblioAnalysis(object(x)))
    t.summary <- e(bibliometrix:::summary.bibliometrix(object = x$overview, k = n, pause = FALSE))
    t.summary <- lapply(t.summary,
      function(x) {
        if (is(x, "data.frame"))
          as_tibble(data.frame(x))
      })
    t.summary <- lst_clear0(t.summary)
    t.summary <- .set_lab(t.summary, sig(x), paste0("data ", names(t.summary)))
    cli::cli_alert_info("bibliometrix:::plot.bibliometrix")
    p.dt <- new_pie(x@tables$step1$alls$DT,
      title = paste0("Total documents: ", nrow(object(x))),
      fun_text = ggrepel::geom_label_repel)
    p.dt <- wrap(p.dt, 10, 8)
    p.dt <- .set_lab(p.dt, sig(x), "Total documents", nrow(object(x)))
    p.most_ave <- .plot_biblioAnalysis(x$overview)
    p.most_ave <- lapply(p.most_ave, function(x) wrap(x, 7, 4))
    p.most_ave <- .set_lab(p.most_ave, sig(x), names(p.most_ave))
    x@plots[[ 2 ]] <- namel(p.most_ave, p.dt)
    x@tables[[ 2 ]] <- namel(t.summary)
    return(x)
  })

setMethod("step3", signature = c(x = "job_biblio"),
  function(x, n = 30, rank.by = "C/Y", journal_info = fxlsx(.prefix("2022-if.xlsx")))
    {
    step_message("Network analysis.")
    require(bibliometrix)
    object(x) <- e(bibliometrix::metaTagExtraction(object(x), Field = "AU_CO"),
      text = "AU_CO: country of author")
    object(x) <- e(bibliometrix::metaTagExtraction(object(x), Field = "AU_UN"),
      text = "AU_UN: affiliation of author")
    if (TRUE) {
      ## country
      obj.net <- e(bibliometrix::biblioNetwork(object(x), analysis = "collaboration",
          network = "countries"))
      p.net <- e(bibliometrix::networkPlot(obj.net, n = dim(obj.net)[1],
          Title = "Country Collaboration", type = "circle", size = TRUE,
          remove.multiple = FALSE, labelsize = 0.8), text = "collaboration")
      p.co_country <- .plot_coNetwork(p.net$graph, fill = "color", size = "deg") +
        guides(fill = "none", edge_width = "none") +
        labs(size = "Centrality degree")
      p.co_country <- wrap(p.co_country, 9, 7)
      p.co_country <- .set_lab(p.co_country, sig(x), "country collaboration network")
    }
    if (TRUE) {
      ## Journal (SO)
      obj <- object(x)
      stat_jour <- head(sort(table(obj$SO), decreasing = TRUE), n = n)
      journal_top <- data.frame(journal = names(stat_jour), count = as.integer(stat_jour))
      match2 <- function(x, y) {
        lst <- lapply(list(x, y),
          function(x) {
            gs(toupper(gs(x, " ", "")), "[^A-Z]", "")
          })
        res <- match(lst[[1]], lst[[2]])
        data.frame(x = x, y = y[res], which = res)
      }
      t.matched_journal_top <- match2(journal_top$journal, journal_info[[1]])
      t.journal_top <- dplyr::slice(journal_info, !!t.matched_journal_top$which)
      t.journal_top <- dplyr::bind_cols(journal_top, t.journal_top)
      t.journal_top <- as_tibble(t.journal_top)
      t.journal_top <- dplyr::select(t.journal_top, journal, count, CATEGORY, `2022IF`)
      t.journal_top <- .set_lab(t.journal_top, sig(x), "top publication journal by count")
    }
    if (TRUE) {
      ## C/Y filter and export (average citation counts per year)
      obj <- object(x)
      obj$C.Y <- obj$TC / (as.integer(format(Sys.time(), "%Y")) - obj$PY + 1)
      if (rank.by == "C/Y") {
        obj <- obj[order(obj$C.Y, decreasing = TRUE), ]
      } else if (rank.by == "TC") {
        obj <- obj[order(obj$TC, decreasing = TRUE), ]
      }
      obj <- head(obj, n = n)
      t.most_citation <- as_tibble(data.frame(obj))
      tags.export <- dplyr::filter(x@params$summary$mandatoryTags,
        status == "Excellent", tag != "NR")
      t.most_citation <- dplyr::select(t.most_citation, short_name = rownames,
        dplyr::all_of(tags.export$tag), citation_per_year = C.Y)
      t.most_citation <- dplyr::rename(t.most_citation,
        dplyr::all_of(nl(tags.export$description, tags.export$tag, FALSE)))
      t.most_citation <- .set_lab(t.most_citation, sig(x), "top ranked documents of citation per year")
    }
    if (TRUE) {
      ## co-citation
      ## by look into the source code of `biblioNetwork`: n, the most frequent reference were selected.
      obj.net <- e(bibliometrix::biblioNetwork(object(x), analysis = "co-citation",
          network = "references"))
      p.net <- e(bibliometrix::networkPlot(obj.net,
          Title = "Co-Citation Network", type = "fruchterman", n = n,
          size = TRUE, remove.multiple = FALSE, labelsize = 0.7, edgesize = 5), text = "co-citation")
      p.co_citation <- .plot_coNetwork(p.net$graph, fill = "color", size = "deg",
        edge_color = "orange", zoRange = 1.5, fun_edge = ggraph::geom_edge_arc,
        width_range = c(.01, .05)
      )
      p.co_citation <- p.co_citation +
        guides(fill = "none", edge_width = "none") +
        labs(size = "Centrality degree")
      p.co_citation <- wrap(p.co_citation, 9, 7)
      p.co_citation <- .set_lab(p.co_citation, sig(x), "top ranked of co citation")
    }
    if (TRUE) {
      ## keywords
      ## all used for keywords analysis
      obj.net <- e(bibliometrix::biblioNetwork(object(x), analysis = "co-occurrences", network = "keywords"))
      ## top n used for plot
      p.net <- e(bibliometrix::networkPlot(obj.net, normalize = "association", weighted = TRUE, n = n,
          Title = "Keyword Co-occurrences", type = "fruchterman", size = TRUE,
          edgesize = 5, labelsize = 0.7), text = "co-occurrences")
      p.co_keywords <- recordPlot()
      p.co_keywords <- .set_lab(p.co_keywords, sig(x), "top ranked keywords of co occurrences")
    }
    x@tables[[ 3 ]] <- namel(t.most_citation, t.journal_top)
    x@plots[[ 3 ]] <- namel(p.co_country, p.co_citation, p.co_keywords)
    return(x)
  })

setMethod("step4", signature = c(x = "job_biblio"),
  function(x, n = 20){
    step_message("Co-Word Analysis and Historical Direct Citation Network.")
    # ID: Keywords Plus associated by ISI or SCOPUS database
    require(bibliometrix)
    # by look into the codes of `bibliometrix::conceptualStructure`,
    # the params `documents` do not results in p.conceptual_structure
    cs <- e(bibliometrix::conceptualStructure(object(x), field = "ID", method = "MCA",
        minDegree = 10, clust = 5, stemming = FALSE, labelsize = 3,
        documents = n, graph = FALSE))
    p.conceptual_structure <- cs$graph_terms +
      theme(
        axis.title = element_text(size = 15, face = "plain"),
        plot.title = element_text(size = 15, face = "plain")
      )
    p.conceptual_structure$layers[[ 6 ]] <- NULL
    p.conceptual_structure <- wrap(p.conceptual_structure, 12, 10)
    p.conceptual_structure <- .set_lab(p.conceptual_structure, sig(x),
      "conceptual structure map of a scientific field")
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
  graph <- create_layout(graph, layout = 'linear', circular = TRUE)
  p <- ggraph(graph) +
    fun_edge(aes(width = !!rlang::sym(width)), color = edge_color, alpha = .5) +
    geom_node_point(aes(size = !!rlang::sym(size), fill = !!rlang::sym(fill)), alpha = .7, shape = 21) +
    geom_node_text(aes(x = x * 1.1,  y = y * 1.1, 
        label = !!rlang::sym(label),
        angle = -((-node_angle(x,  y) + 90) %% 180) + 90), 
      size = 3, hjust = 'outward', family = "Times") +
    scale_size(range = c(3, 9)) +
    scale_fill_manual(values = color_set()) +
    scale_edge_width_continuous(range = width_range) +
    coord_cartesian(xlim = zoRange(graph$x, zoRange), ylim = zoRange(graph$y, zoRange)) +
    theme_graph() +
    theme(text = element_text(family = "Times"))
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

