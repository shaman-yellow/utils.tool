# ==========================================================================
# workflow of pubmed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_pubmed <- setClass("job_pubmed", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "pubmed",
    info = c("https://www.ncbi.nlm.nih.gov/books/NBK3837/"),
    cite = "",
    method = "",
    tag = "pubmed",
    analysis = "PubMed 数据搜索"
    ))

job_pubmed <- function(keys, type = c("Review", "Clinical Study"),
  elements = c("SortPubDate", "Title", "FullJournalName", "Id"),
  ...)
{
  query <- add_field(keys, "[Title/Abstract]")
  type <- match.arg(type)
  extra <- glue::glue(
    "{type}[Filter]"
  )
  object <- edirect_db("pubmed", c(query, extra), elements, ...)
  object <- .set_lab(object, keys[1], "EDirect query")
  x <- .job_pubmed(object = object, params = namel(query, elements))
  x <- methodAdd(x, "使用 Entrez Direct (EDirect) <https://www.ncbi.nlm.nih.gov/books/NBK3837/> 搜索 PubMed 数据库 (`esearch -db pubmed`)，查询信息为: {query}。")
  x
}

setMethod("step0", signature = c(x = "job_pubmed"),
  function(x){
    step_message("Prepare your data with function `job_pubmed`.")
  })

pubmed_types <- function() {
  c("Adaptive Clinical Trial",
    "Address",
    "Autobiography",
    "Bibliography",
    "Biography",
    "Books and Documents",
    "Case Reports",
    "Classical Article",
    "Clinical Conference",
    "Clinical Study",
    "Clinical Trial",
    "Clinical Trial Protocol",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Clinical Trial, Veterinary",
    "Collected Work",
    "Comment",
    "Comparative Study",
    "Congress",
    "Consensus Development Conference",
    "Consensus Development Conference, NIH",
    "Controlled Clinical Trial",
    "Corrected and Republished Article",
    "Dataset",
    "Dictionary",
    "Directory",
    "Duplicate Publication",
    "Editorial",
    "Electronic Supplementary Materials",
    "English Abstract",
    "Equivalence Trial",
    "Evaluation Study",
    "Expression of Concern",
    "Festschrift",
    "Government Publication",
    "Guideline",
    "Historical Article",
    "Interactive Tutorial",
    "Interview",
    "Introductory Journal Article",
    "Lecture",
    "Legal Case",
    "Legislation",
    "Letter",
    "Meta-Analysis",
    "Multicenter Study",
    "News",
    "Newspaper Article",
    "Observational Study",
    "Observational Study, Veterinary",
    "Overall",
    "Patient Education Handout",
    "Periodical Index",
    "Personal Narrative",
    "Portrait",
    "Practice Guideline",
    "Pragmatic Clinical Trial",
    "Preprint",
    "Published Erratum",
    "Randomized Controlled Trial",
    "Randomized Controlled Trial, Veterinary",
    "Research Support, American Recovery and Reinvestment Act",
    "Research Support, N.I.H., Extramural",
    "Research Support, N.I.H., Intramural",
    "Research Support, Non-U.S. Gov't",
    "Research Support, U.S. Gov't, Non-P.H.S.",
    "Research Support, U.S. Gov't, P.H.S.",
    "Research Support, U.S. Gov't",
    "Retracted Publication",
    "Retraction of Publication",
    "Review",
    "Scientific Integrity Review",
    "Systematic Review",
    "Technical Report",
    "Twin Study",
    "Validation Study",
    "Video-Audio Media",
    "Webcast"
  )
}
