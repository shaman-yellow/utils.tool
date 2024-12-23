# ==========================================================================
# Fast display the content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

send_eval <- function(to,
  subject = "需要评估绩效的业务",
  content = "Hello, 慧姐\n\n这是这个月需要评估绩效的业务 (各个业务 ICDTM 的汇总，如有系数为0，则为空置，待评估)，已附在表格中。\n\nBest wishes!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0(.prefix(), "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "need_eval.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_registers <- function(to,
  subject = "业务表格更新",
  content = "Hello, 慧姐\n\n这是每月末需提交的更新的业务登记表。\n\nBest wishes!",
  path_summary = .prefix("summary"),
  atts = c(file.path(path_summary, "生信组表格登记-黄礼闯.xlsx"),
    file.path(path_summary, "type_list.xlsx")))
{
  send_that(to, subject, content, atts)
}

send_prin <- function(to,
  subject = "月底交个人行为准则考核表",
  content = "Hello, 慧姐\n\n个人行为准则考核表已交。\n\nBest wishes!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0(.prefix(), "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "2023年行为准则考核表.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_summary <- function(to,
  subject = "月底交付考核",
  content = "Hello, 慧姐\n\n绩效表已提交。\n\nBest wishes!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0(.prefix(), "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "assess_绩效+软性考核表.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_that <- function(to, subject, content, atts = NULL, from = "huanglichuang@wie-biotech.com") {
  isthat <- usethis::ui_yeah("Are you sure sending the mail?")
  if (!isthat) {
    return(message("Sending cancelled."))
  }
  if (is.character(content)) {
    file <- tempfile("mail_", fileext = ".txt")
    writeLines(content, file)
  }
  if (!is.null(atts)) {
    code_atts <- paste0(" -a ", paste0(atts, collapse = " "), " -- ")
  } else {
    code_atts <- ""
  }
  cdRun("mutt ", to,
    " -e 'my_hdr From:", from, "'",
    " -s '", subject, "'",
    code_atts,
    " < ", file
  )
  message("Done")
}

update_tagsTable <- function(tags = .tag_anno(),
  path = .prefix("summary"), file = "type_list.xlsx")
{
  tags.cn <- unique(tags)
  data <- tibble::tibble(Seq = seq_along(tags.cn), Name = tags.cn)
  data <- tibble::add_row(data,
    Name = c("...",
      paste0("Last update: ", Sys.time()),
      "Author: Huang LiChuang")
  )
  file <- paste0(path, "/", file)
  openxlsx::write.xlsx(data, file)
  browseURL(normalizePath(file))
}

update_registers <- function(orders = get_orders(),
  path = .prefix("summary"),
  name = "黄礼闯",
  register = paste0(path, "/", "生信组表格登记-", name, ".xlsx"),
  templ_dir = .prefix("performance_templ/"))
{
  target <- paste0(templ_dir, "/", "生信组表格登记-XX.xlsx")
  if (!dir.exists(path)) {
    dir.create(path)
  }
  if (!file.exists(register)) {
    file.copy(target, register)
  }
  wb <- openxlsx2::wb_load(register)
  pos <- list(2, 1)
  ## base
  data <- dplyr::filter(orders, type != "备单业务")
  data <- dplyr::mutate(data, seq = seq_len(nrow(data)), note = as.note.tags(tags.cn))
  data <- dplyr::select(data,
    date, seq, info, id, score, member, receive_date, status, title, note
  )
  wb <- openxlsx2::wb_add_data(wb, 1, data, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos))
  ## extra
  data <- dplyr::filter(orders, type == "备单业务")
  data <- dplyr::mutate(data, seq = seq_len(nrow(data)), note = as.note.tags(tags.cn))
  data <- dplyr::select(data,
    date, seq, info, id, score, member, receive_date, status, title, note
  )
  wb <- openxlsx2::wb_add_data(wb, 2, data, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos))
  openxlsx2::wb_save(wb, register)
  browseURL(normalizePath(register))
}

summary_month <- function(
  time = Sys.time(),
  orders = get_orders(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path = .prefix(),
  templ_dir = .prefix("performance_templ/"),
  rm = F)
{
  dir <- paste0(path, "/", "summary_", paste0(year, "-", month))
  targets <- c(ass = "assess_绩效+软性考核表.xlsx", prin = "2023年行为准则考核表.xlsx")
  if (rm) {
    unlink(list.files(dir, full.names = T, all.files = T, recursive = T), T, T)
    dir.create(dir)
    file.copy(paste0(templ_dir, "/", targets), dir)
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
    file.copy(paste0(templ_dir, "/", targets), dir)
  }
  targets[] <- paste0(dir, "/", targets)
  wb <- openxlsx2::wb_load(targets[[ "ass" ]])
  ## prepare orders (this month)
  orders <- dplyr::filter(orders, lubridate::year(belong) == !!year,
    lubridate::month(belong) == !!month)
  ## modify the assess table
  pos.date <- list(3, 1)
  date <- paste0(year, "年 ", month, "月 01日", "至",
    year, "年 ", month, "月 ", lubridate::days_in_month(time), "日")
  wb <- openxlsx2::wb_add_data(wb, 1, date,
    dims = do.call(openxlsx2::wb_dims, pos.date))
  if (T) {
    pos.data_ass <- list(7, 1)
    data_ass <- dplyr::mutate(orders,
      seq = seq_len(nrow(orders)), num = 1, title.en = "",
      note = as.note.tags(tags.cn),
      coef = round(coef, 3)
    )
    data_ass <- dplyr::select(data_ass,
      member, seq, id, type, class, score, num, title, title.en, status, note, coef,
      remuneration, I, C, D, `T`, M
    )
    fun <- function(wb, data) {
      openxlsx2::wb_add_data(wb, 1, data, col_names = F,
        dims = do.call(openxlsx2::wb_dims, pos.data_ass), na.strings = "")
    }
    CopyNeedEval <- F
    if (any(is.na(data_ass$coef))) {
      wb_eval <- fun(wb, dplyr::filter(data_ass, is.na(coef)))
      openxlsx2::wb_save(wb_eval, paste0(dir, "/need_eval.xlsx"))
    } else {
      CopyNeedEval <- T
    }
    wb <- fun(wb, data_ass)
  }
  if (T) {
    pos.data_sum <- list(29, 6)
    data_ass <- dplyr::mutate(data_ass,
      type = ifelse(type %in% c("固定业务", "备单业务"), "base", "other"))
    if (!any(data_ass$type == "base")) {
      data_ass <- tibble::add_row(data_ass, type = "base")
    }
    data_ass <- dplyr::arrange(data_ass, type)
    data_sum <- split_lapply_rbind(data_ass, ~type,
      function(data) {
        if (any(is.na(data$coef))) {
          return(data.frame(num = 0, score = "", coef = 0, coef_fomu = ""))
        }
        sum <- data.frame(
          num = nrow(data),
          score = "",
          coef = sum(data$coef)
        )
        dplyr::mutate(sum,
          coef_fomu = paste0(paste0(data$coef, collapse = "+"), "=", coef)
        )
      })
    wb <- openxlsx2::wb_add_data(wb, 1, data_sum, col_names = F,
      dims = do.call(openxlsx2::wb_dims, pos.data_sum))
  }
  if (T) {
    pos.sum <- list(29, 10)
    sum <- sum(data_sum$coef)
    wb <- openxlsx2::wb_add_data(wb, 1, sum,
      dims = do.call(openxlsx2::wb_dims, pos.sum))
  }
  openxlsx2::wb_save(wb, targets[[ "ass" ]])
  if (CopyNeedEval) {
    file.copy(targets[[ "ass" ]], paste0(dir, "/need_eval.xlsx"))
  }
  ########################################
  ########################################
  wb <- openxlsx2::wb_load(targets[[ "prin" ]])
  info <- data.frame(depart = "部门：技术部",
    x1 = "", x2 = "",
    job = "岗位：生物信息工程师",
    x3 = "", name = "姓名：黄礼闯", x4 = "",
    month = paste0("月度：", year, "年", month, "月"))
  pos.info <- list(2, 1)
  wb <- openxlsx2::wb_add_data(wb, 1, info, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos.info))
  pos.sign <- list(22, 1)
  sign <- paste0("本人确认：            日 期 ：", year, "年 ", month, "月  日")
  wb <- openxlsx2::wb_add_data(wb, 1, sign,
    dims = do.call(openxlsx2::wb_dims, pos.sign))
  openxlsx2::wb_save(wb, targets[[ "prin" ]])
  ## check
  browseURL(normalizePath(targets[[ "ass" ]]))
}

.thedir_job <- function(month, year = 2023, path = .prefix("backup", "db")) {
  month <- as.character(month)
  if (nchar(month) < 2) {
    month <- paste0("0", month)
  }
  paste0(path, "/", year, month)
}

.cp_job <- function(from, to) {
  if (!dir.exists(to)) {
    dir.create(to)
  }
  file.copy(from, to, T, T)
}

backup_jobs <- function(alls, time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time))
{
  data <- dplyr::filter(alls,
    lubridate::month(belong) == !!month,
    lubridate::year(belong) == !!year
  )
  thedir <- .thedir_job(month, year)
  num <- 0L
  dir.create(thedir, F)
  res <- lapply(data$.dir,
    function(dir) {
      files <- list.files(dir, ".*\\.zip", full.names = T)
      if (length(files) > 1) {
        n <- menu(basename(files), title = "Select a zip file to backup.")
        files <- files[n]
      }
      if (!length(files)) {
        message("No zip files found in dir: ", dir)
        isThat <- usethis::ui_yeah("Select other files?")
        if (isThat) {
          files <- list.files(dir, full.names = T)
          files <- files[ menu(files, title = "All available files") ]
        }
      }
      if (length(files)) {
        num <<- num + 1L
        .cp_job(files, paste0(thedir, "/", num, "_", basename(files)))
      }
    })
  unlist(res)
}

move_rds_all <- function(alls, time, to = "~/disk_sdb1")
{
  if (!dir.exists(to)) {
    stop("dir.exists(to)")
  }
  data <- dplyr::filter(alls, belong %in% as.Date(time))
  lapply(data$.dir, function(x) move_rds(x, to))
}

move_raw_all <- function(alls, time, before = F, to = "~/disk_sdb1",
  prefix = c("GDC", "GSE", "qiime", "sra"))
{
  if (before)
    data <- dplyr::filter(alls, belong <= as.Date(time))
  else
    data <- dplyr::filter(alls, belong %in% as.Date(time))
  pbapply::pblapply(data$.dir,
    function(dir) {
      path <- grpf(list.dirs(dir), paste0(paste0(prefix, "[^/]+$"), collapse = "|"))
      path <- path[ !dirname(path) %in% path ]
      if (length(path)) {
        target <- paste0(to, "/", basename(gs(dir, "/$", "")))
        if (!dir.exists(target)) {
          dir.create(target)
        }
        res <- file.copy(path, target, T, T)
        if (all(unlist(res)))
          unlink(path, recursive = T, force = T)
      }
    })
}

move_rds <- function(from, to, exclude = c(".items.rds", ".injobs.rds")) {
  lapply(from,
    function(path) {
      dirname <- basename(gs(path, "/$", ""))
      target <- paste0(to, "/", dirname)
      dir.create(target, F)
      files <- list.files(path, pattern = ".rds$|.RData$|.Rdata$|.rdata$", full.names = T, all.files = T)
      thefilenames <- basename(files)
      files <- files[ !thefilenames %in% exclude ]
      lapply(files,
        function(file) {
          if (.Platform$OS.type == "unix") {
            system(paste0("cp ", file, " -t ", target))
            file.remove(file)
          }
        })
    })
}

cf <- function(remuneration, base_wage = 6000) {
  remuneration / base_wage
}

odate <- function(month, year = format(Sys.time(), "%Y")) {
  as.Date(paste0(year, "-", month, "-01"))
}

.get_wages <- function(date, base_wage) {
  vapply(date, FUN.VALUE = double(1),
    function(date) {
      which <- utils::tail(which(as.Date(names(base_wage)) <= as.Date(date)), n = 1)
      base_wage[[which]]
    })
}

.remu <- function(coef, date, base_wage) {
  use <- .get_wages(date, base_wage)
  coef * use
}

.setFill.orders <- function() {
  list("2023-12-01" = 3800)
}

thismonth <- function(data = get_orders(),
  month = lubridate::month(time), year = lubridate::year(time),
  time = Sys.Date())
{
  month <- as.Date(paste0(year, "-", month, "-01"))
  data <- dplyr::filter(data, belong == month)
  writeLines(paste0("\nCoef: ", sum(data$coef), "\nRem: ", sum(data$remuneration)))
  invisible(data)
}

get_orders <- function(
  dir = .prefix(), pattern = "\\.items.*\\.rds",
  base_wage = list("2023-07-01" = 3000, "2023-08-01" = 4500, "2023-09-01" = 6000),
  setFill = .setFill.orders())
{
  files <- list.files(dir, pattern, T, T, T)
  lst <- lapply(files,
    function(file) {
      lst <- readRDS(file)
      maybeMulti <- c("coef", "belong", "id", "type", "eval", "I", "C", "D", "T", "M")
      maybeMulti <- names(lst[ names(lst) %in% maybeMulti & lengths(lst) > 1 ])
      lst.m <- lst[ names(lst) %in% maybeMulti ]
      lst <- lst[ !names(lst) %in% maybeMulti ]
      lst <- lapply(lst, function(x) if (identical(x, character(0))) character(1) else x)
      data <- do.call(tibble::tibble, lst)
      if (length(lst.m$coef) > 1) {
        data <- do.call(rbind, rep(list(data), length(lst.m$coef)))
      }
      data <- do.call(dplyr::mutate, c(list(data), lapply(lst.m, unlist)))
      data <- dplyr::mutate(data,
        coef = as.double(coef),
        date = as.Date(date),
        receive_date = as.Date(receive_date),
        .dir = dirname(file),
        tags.list = lapply(strsplit(tags, ","), gs, "\\s", ""),
        tags.cn = recode_tags(tags.list)
      )
      data
    })
  data <- as_tibble(data.table::rbindlist(lst, fill = T))
  data <- dplyr::arrange(data, belong, receive_date)
  data <- dplyr::mutate_if(data, is.character,
    function(x) ifelse(is.na(x), "", x))
  if (!is.null(setFill)) {
    for (i in seq_along(setFill)) {
      i.month <- names(setFill)[i]
      if (is.na(as.Date(i.month, optional = T))) {
        i.month <- paste0(i.month, "-01")
      }
      theWhich <- which(grpl(data$belong, paste0("^", i.month)) & is.na(data$coef))
      if (length(theWhich)) {
        i.base <- .get_wages(i.month, base_wage)
        data$coef[theWhich] <- cf(setFill[[ i ]] / length(theWhich), i.base)
      }
    }
  }
  data <- dplyr::mutate(data, remuneration = .remu(coef, belong, base_wage))
  data <- dplyr::relocate(data, belong, id, title, remuneration)
  print(data, n = Inf)
  invisible(data)
}

recode_tags <- function(lst) {
  fun <- function(tags) {
    dic <- as.list(.tag_anno())
    tags <- vapply(tags, function(x) if (identical(x, character(0))) "" else x, character(1))
    notIn <- (!tags %in% names(dic)) & (tags != "")
    if (any(notIn)) {
      message("Function: .tag_anno()")
      stop("Some tags without annotation:\n\t", paste0(tags[ notIn ], collapse = ", "))
    }
    res <- dplyr::recode(tags, !!!dic, .default = character(1))
    res
  }
  if (is(lst, "list")) {
    lapply(lst, fun)
  } else {
    fun(lst)
  }
}

as.note.tags <- function(tags) {
  fun <- function(x) {
    x <- x[ x != "" ]
    if (length(x)) {
      paste0(paste0("(", seq_along(x), ") ", x), collapse = "; ")
    } else {
      ""
    }
  }
  if (is(tags, "list")) {
    vapply(tags, fun, character(1))
  } else {
    fun(tags)
  }
}

plot_orders_summary <- function(data) {
  data <- dplyr::mutate(data,
    .type = dplyr::recode(type, `备单业务` = "BI", `固定业务` = "IN", "其他业务" = "Others")
  )
  message("Plot all proportion")
  p.all.pop <- new_pie(data$.type, title = "All types")
  idata <- split(data, ~belong)
  fun <- function(lst) {
    lst <- lapply(lst,
      function(data) {
        dplyr::mutate(data, value = 1 / nrow(data))
      })
    data <- do.call(dplyr::bind_rows, lst)
    p <- ggplot(data, aes(x = reorder(gs(belong, "-01$", ""), belong, decreasing = TRUE), y = value)) +
      geom_col(aes(fill = .type), position = "stack", width = .6) +
      coord_flip() +
      labs(y = "Months", x = "Proportion", fill = "Type") +
      scale_fill_manual(values = ggsci::pal_npg()(8)) +
      theme_minimal()
    p
  }
  message("Plot individual proportion")
  p.ind.pop <- fun(idata)
  fun <- function(data) {
    data <- dplyr::summarize(dplyr::group_by(data, belong), sum_coef = sum(coef))
    p <- ggplot(data, aes(x = seq_len(nrow(data)), y = sum_coef)) +
      geom_line() +
      geom_area(fill = 'grey90', alpha = .5) +
      scale_x_continuous(labels = gs(data$belong, "-01$", ""),
        breaks = unique(as.integer(data$belong))) +
      coord_cartesian(ylim = zoRange(data$sum_coef, 5)) +
      labs(x = "Month", y = "Coef") +
      theme_minimal()
    wrap(p, 5, 3)
  }
  message("Plot coefficients")
  p.all.coef <- fun(data)
  fun <- function(data) {
    data <- dplyr::mutate(data, used = date - receive_date + 1)
    data <- dplyr::summarize(dplyr::group_by(data, .type, belong), avg_time = mean(used))
    data <- dplyr::arrange(data, belong)
    p <- ggplot(data, aes(x = as.integer(belong), y = avg_time, color = .type)) +
      geom_line() +
      ggsci::scale_color_npg() +
      scale_x_continuous(labels = gs(unique(data$belong), "-01$", ""),
        breaks = unique(as.integer(data$belong))) +
      labs(x = "Month", y = "Cost time (Day)", color = "Type") +
      theme_minimal()
    wrap(p, 7, 3)
  }
  message("Plot time used.")
  p.ind.avgTime <- fun(data)
  null <- nullGrob()
  fr1 <- frame_col(c(p.all.pop = 1, p.ind.pop = 1),
    list(p.all.pop = p.all.pop@data, p.ind.pop = as_grob(p.ind.pop)))
  fr2 <- frame_col(c(p.all.coef = 1, p.ind.avgTime = 1),
    list(p.all.coef = as_grob(p.all.coef@data),
      p.ind.avgTime = as_grob(p.ind.avgTime@data)))
  p.alls <- wrap(yf(fr1 = 1, fr2 = 1))
  p.alls <- wrap(p.alls, 10, 5)
  namel(p.all.pop, p.ind.pop, p.all.coef, p.ind.avgTime, p.alls)
}

plot_task_summary <- function() {
  route <- as_network(list(
      "Mail:Deparse",
      "Deparse:Tags",
      "Tags:Info, Analysis",
      "Info:Extracting",
      "Extracting:ID, Score, Client, Date..",
      "ID, Score, Client, Date..:Namespace",
      "Analysis:Workflow",
      "Workflow:Job.1, job.2, job...",
      "Job.1, job.2, job...:Summary",
      "Summary:Report",
      "Report:Packaging",
      "Packaging, Namespace:Record",
      "Packaging:Send_mail"
      ))
  p.each <- flowChart(route, 1.1, 1)
  p.each <- wrap(p.each, 7.5, 7)
  route <- as_network(list(
      "Monthly_Tasks:Task.1, Task.2, Task.3, Task...",
      "Task.1, Task.2, Task.3, Task...:Summary",
      "Summary:Titles, Scores, IDs, Dates, Analysis",
      "Titles, Scores, IDs, Dates, Analysis:Information",
      "Information:Tables",
      "Tables:Record_1, Record_2, Record_3, Record_4",
      "Record_1, Record_2, Record_3, Record_4:Send_mail"
      ))
  p.month <- flowChart(route, 1.1, 1)
  p.month <- wrap(p.month, 7.5, 7)
  p.alls <- frame_col(c(p.each = 1.3, p.month = 1),
    list(p.each = as_grob(p.each@data),
      p.month = as_grob(p.month@data)))
  p.alls <- wrap(p.alls, 12, 7)
  namel(p.each, p.month, p.alls)
}

plot_report_summary <- function(cover_file = .prefix("cover_page.pdf"))
{
  outline <- function(name) {
    grecti(name, tfill = "black",
      tgp_args = list(col = "transparent"), borderF = 1.5)
  }
  block <- function(text) {
    into(grectn("grey92", , list(col = "transparent")),
      gtext(text, list(font = c(plain = 1)), x = .1, hjust = 0))
  }
  stext <- function(text) {
    gtext(text, x = 0, hjust = 0)
  }
  txt.h1 <- stext("Table name")
  txt.h2 <- stext("Table location")
  blk.b1 <- block("Dim: w * h")
  blk.b2 <- block("- Col1: ...\n- Col2: ...")
  df <- data.frame(Col1 = n(id., 3), Col2 = n(cand., 3))
  grob_table <- gridExtra::tableGrob(df, theme = gridExtra::ttheme_default(8))
  null <- nullGrob()
  tab.cont <- yf(txt.h1 = 1, txt.h2 = 1, null = .1, blk.b1 = 1, null = .2,
    blk.b2 = 2, null = .2, grob_table = 4)
  tab.cont <- ggather(tab.cont, vp = viewport(, , .8, .8))
  tab <- outline("Table show")
  tab <- into(tab, tab.cont)
  ## figure
  fig <- into(grectn(), gtext("Fig.", list(cex = 3)))
  txt.h1 <- stext("Table name")
  txt.h2 <- stext("Table location")
  fig.cont <- yf(null = 1, txt.h1 = 1, txt.h2 = 1, null = .2, fig = 3, null = 1.5)
  fig.cont <- ggather(fig.cont, vp = viewport(, , .8, .8))
  fig <- outline("Table show")
  fig <- into(fig, fig.cont)
  ## table of contents
  fun <- function(text, name) {
    obj <- lapply(1:3, function(x) stext(paste0("\t", text, " ", x)))
    names(obj) <- paste0(name, 1:3)
    obj
  }
  ans <- fun("Analysis", "ans")
  ts <- fun("Table", "ts")
  fs <- fun("Figure", "fs")
  txt.t1 <- stext("Contents")
  txt.t2 <- stext("Tables")
  txt.t3 <- stext("Figures")
  wei.lst <- nl(c("txt.t1", names(ans), "txt.t2", names(ts), 'txt.t3', names(fs)),
    1, F, default = 1)
  lst.cont <- frame_row(wei.lst, c(ans, ts, fs, namel(txt.t1, txt.t2, txt.t3, null)))
  lst.cont <- ggather(lst.cont, vp = viewport(, , .8, .8))
  lst <- outline("Table of contents")
  lst <- into(lst, lst.cont)
  ## Cover
  tmp <- tempfile()
  pdf_convert(cover_file, filenames = tmp)
  cover <- ggather(rasterGrob(png::readPNG(tmp)), vp = viewport(, , .8, .9))
  cover <- into(outline("Cover"), cover)
  p.rep <- xf(cover = 1, null = .1, lst = 1, null = .1, tab = 1, null = .1, fig = 1)
  wrap(p.rep, 9, 5)
}

generate_tiffs <- function(pattern = "^MAIN-Fig.*\\.pdf",
  dir = get_savedir("figs"), output = "TIFF")
{
  output <- file.path(dir, output)
  dir.create(output, F)
  files <- list.files(dir, pattern, full.names = T)
  lapply(files,
    function(x) {
      newfile <- file.path(output, paste0(tools::file_path_sans_ext(basename(x)), ".tiff"))
      pdf_convert(x, "tiff", filenames = newfile, dpi = 300)
    })
  return(output)
}

order_packaging <- function(target = "output.pdf",
  register = autoRegisters, idname = gidn(), ...)
{
  report <- paste0(idname, ".", tools::file_ext(target))
  if (!file.exists(report)) {
    file.copy(target, report, T)
  }
  if (!is.null(relocate <- getOption("autoRegisters_relocate"))) {
    if (!identical(names(register), names(relocate))) {
      message('identical(names(register), names(relocate)) == F, not match in names.')
      register <- register[ names(register) %in% names(relocate) ]
    }
    lapply(names(register), function(name) {
      dir <- dirname(relocate[[ name ]])
      dir.create(dir, F, recursive = T)
      file.copy(register[[name]], relocate[[name]], T, T)
    })
    tempfiles <- unlist(as.list(relocate))
    package_results(main = tempfiles, head = NULL, masterZip = NULL, report = report, ...)
    unlink(tempfiles, T, T)
  } else {
    package_results(main = register, head = NULL, masterZip = NULL, report = report, ...)
  }
  zipfiles <- list.files(".", "[\u4e00-\u9fa5]+.*\\.zip")
  zipfiles <- zipfiles[ zipfiles != paste0(idname, ".zip") ]
  if (length(zipfiles)) {
    isThat <- usethis::ui_yeah(
      paste0("Detected exists zip files:\n\t", paste0(zipfiles, "\n\t"), "\nRemove thats?")
    )
    if (isThat) {
      file.remove(zipfiles)
    }
  }
  file.rename("./client.zip", paste0(idname, ".zip"))
  browseURL(paste0(idname, ".zip"))
  return(idname)
}

package_results <- function(
  master_include = NULL,
  report = "output.pdf",
  main = unique(c(search_resFile("index.Rmd"), autoRegisters)), 
  headPattern = "[^r][^e][^s].*业务.*\\.zip", head = list.files(".", headPattern),
  masterSource = c("output.Rmd", "install.R", "read_me.txt", head),
  prefix = "results_", zip = paste0(prefix, head),
  clientZip = "client.zip", masterZip = "master.zip",
  clear = T, external_file = "order_material", extras = NULL)
{
  if (!is.null(external_file)) {
    if (!file.exists(external_file)) {
      stop("file.exists(external_file)")
    } else {
      if (length(exters <- list.files(external_file, full.names = T)) == 0) {
        stop("Guess you have forget to save the order material file to `external_file`.")
      } else {
        dir.create(extern0 <- "报单相关资料", F)
        file.copy(exters, extern0, T, T)
      }
    }
  } else {
    extern0 <- NULL
  }
  files.client <- c(report, main, extern0, extras)
  files.client <- files.client[ !duplicated(normalizePath(files.client)) ]
  files.master <- c(files.client, head, masterSource)
  zip(clientZip, files.client)
  if (!is.null(masterZip)) {
    zip(masterZip, files.master)
    zip(zip, c(clientZip, masterZip), flags = "-ur9X")
    if (clear) {
      file.remove(c(clientZip, masterZip))
    }
  } else {
    message("Only client file were packaged.")
  }
  if (file.exists(file <- ".items.rds")) {
    info <- readRDS(file)
    info$date <- Sys.Date()
    info$isLatest <- T
    saveRDS(info, file)
  }
}

od_get_title <- function() {
  odb("name", "analysis")
}

cal_coef <- function(x, a = .03, k = .028) {
  if (is(x, "ic")) {
    x <- sum(unlist(x))
  } else if (is(x, "list")) {
    x <- vapply(x, function(x) sum(unlist(x)), double(1))
  }
  a + x * k
}

ic <- function(Inexperience = 0, Complexity = 0, Deficiency = 0, Time_consuming = 0, Massiveness = 0)
{
  .ic(list(I = Inexperience, C = Complexity, D = Deficiency, 'T' = Time_consuming, M = Massiveness))
}

show.ic <- function(info, a = .03, k = .028, base_wage = 6000)
{
  data <- data.frame(I = info$I, C = info$C, D = info$D, 'T' = info$T, M = info$M)
  print(knitr::kable(data))
  str <- paste0("\n", "y = a + k * (I + C + D + T + M) = ",
    a, " + ", k, " * (",
    info$I, " + ", info$C, " + ", info$D, " + ", info$T, " + ", info$M,
    ") = ", info$coef, "", "\n"
  )
  writeLines(str)
}

.ic <- setClass("ic", 
  contains = c("list"),
  representation = representation(),
  prototype = NULL)

get_all_tags <- function(rm = NULL, jobs = .get_job_list()) {
  if (length(jobs)) {
    tags <- vapply(jobs, FUN.VALUE = character(1),
      function(x) {
        tags <- ""
        if (is(x, "character")) {
          x <- get(x, envir = .GlobalEnv)
        }
        if (is(x, "job")) {
          if (!is.null(rm)) {
            if (class(x) %in% rm) {
              return("")
            }
          }
          if (inherits(try(validObject(x), T), "try-error")) {
            x <- upd(x)
          }
          tags <- x@tag
        }
        tags
      })
    paste0(unique(unlist(strsplit(tags, ", "))), collapse = ", ")
  } else {
    return("")
  }
}

items <- function(
  belong = as.Date(receive_date),
  coef = cal_coef(eval),
  eval = ic(),
  type = od_guess_type(),
  title = od_get_title(),
  status = NA_character_,
  date = Sys.Date(),
  start = NA,
  end = NA,
  finish = NA,
  client = "",
  inst = "",
  info = od_get_info(),
  id = od_get_id(),
  receive_date = od_get_date(),
  score = od_get_score(),
  member = "黄礼闯",
  save = ".items.rds",
  tags = get_all_tags(),
  note = character(1),
  class = "生信分析",
  isLatest = F,
  lock = F)
{
  if (length(class) > 1) {
    if (is.null(getOption("orderClass"))) {
      message("Class ...")
      class <- class[menu(class, title = "Is which type?")]
      options(orderClass = class)
    } else {
      class <- getOption("orderClass")
    }
  }
  if (missing(type)) {
    if (!missing(coef)) {
      type <- ifelse(vapply(coef, function(x) identical(x, .25), logical(1)), "固定业务", "其他业务")
    }
  }
  if (!missing(coef)) {
    if (any(is.na(coef))) {
      which <- which(is.na(coef))
      if (is(eval, "list")) {
        if (is(eval[[ which ]], "ic")) {
          coef[which] <- cal_coef(eval[[ which ]])
          rm(which)
        }
      }
    }
  }
  if (is.null(id)) {
    stop("The `id` can not be a NULL !!!")
  }
  if (identical(id, "")) {
    stop("The `id` can not be a empty character !!!")
  }
  if (is.null(title)) {
    stop("is.null(title)")
  }
  items <- as.list(environment())
  if (T) {
    # new feature
    eval <- items$eval
    items <- items[ names(items) != "eval" ]
    items$icEval <- T
    if (is(eval, "ic")) {
      items <- c(items, as.list(eval))
    } else if (is(eval, "list")) {
      n <- 0L
      eval <- lapply(eval[[1]],
        function(x) {
          n <<- n + 1
          vapply(seq_along(eval), function(w) eval[[ w ]][[ n ]], double(1))
        })
      items <- c(items, as.list(eval))
    }
  }
  if (file.exists(save)) {
    if (!usethis::ui_yeah("File exists, use this?")) {
      return(writeLines("Stopped."))
    }
    info <- readRDS(save)
    if (is.null(info$isLatest)) {
      items$date <- info$date
    } else if (info$isLatest) {
      items$date <- info$date
    }
  }
  if (is.null(lock) || !lock) {
    saveRDS(items, save)
  } else {
    warning("This object has been manually modified, can not update automatically.")
  }
  items
}

od_guess_type <- function() {
  lines <- od_get(pattern = NULL)
  if (any(grpl(lines, "固定业务"))) {
    "固定业务"
  } else {
    "其他业务"
  }
}

gidn <- gid <- function(theme = NULL, items = info, member = 3) {
  idn <- character(0)
  if (grpl(items$id, "^[A-Za-z]*[0-9]+$")) {
    idn <- paste0(items$id, "-", member)
  } else {
    id <- odk("id")
    if (!is.null(id)) {
      idn <- id
    }
  }
  ext <- function(key) {
    stringr::str_extract(items$info, paste0("(?<=", key, "：)", "[^ ]+"))
  }
  sale <- ext("销售")
  if (!length(sale)) {
    sale <- NA
  }
  if (is.na(sale)) {
    sale <- odk("sale")
  }
  if (!is.null(sale)) {
    idn <- paste0(idn, "+销售：", sale)
  }
  client <- ext("客户")
  if (!length(client)) {
    client <- NA
  }
  if (is.na(client)) {
    client <- odk("client")
    if (is.null(client)) {
      client <- .guess_client(od_get(pattern = NULL))
    }
  }
  if (!is.null(client)) {
    if (identical(idn, character(1)) || identical(idn, character(0))) {
      if (identical(client, character(1)) || identical(client, character(0))) {
        stop("Client not setup.")
      }
      idn <- paste0(client, "订单+", odk("analysis"))
    } else {
      idn <- paste0(idn, "+客户：", client)
    }
  } else {
    if (identical(idn, character(0))) {
      idn <- odb("client", "analysis", collapse = "+")
    }
    if (idn == "") {
      stop("Cant not found the mark of `client` or `analysis`")
    }
  }
  if (!is.null(theme)) {
    idn <- paste0(idn, "+", theme)
  } else {
    theme <- odk("theme")
    if (!is.null(theme)) {
      idn <- paste0(idn, "+", theme)
    } else {
      theme <- odk("analysis")
      if (is.null(theme)) {
        stop("No 'analysis[[...]] tags found.'")
      }
      if (!grpl(idn, theme))
        idn <- paste0(idn, "+", theme)
    }
  }
  if (!is.na(items$score)) {
    if (items$score != "") {
      if (grpl(items$score, "[0-9]")) {
        idn <- paste0(idn, "+", items$score, "分")
      }
    }
  }
  idn
}

od_get_id <- function(reload = F, ...) {
  if (reload)
    ID <- NULL
  else
    ID <- getOption("order_id")
  if (is.null(ID)) {
    it <- od_get(..., key = "id")
    if (is.null(it) || it == "") {
      it <- odk("id")
      if (is.null(it))
        it <- odb("client", "analysis")
    }
    options(order_id = it)
    it
  } else {
    ID
  }
}

od_get_score <- function(...) {
  it <- od_get(..., key = "score")
  if (is.null(it) || it == "") {
    lines <- od_get(..., key = "info", pattern = NULL)
    ext <- function(key) {
      res <- unlist(stringr::str_extract_all(lines, paste0("(?<=", key, "：)", ".+")),
        use.names = F)
      res <- res[!is.na(res)]
      res[ nchar(res) == max(nchar(res)) ][1]
    }
    it <- ext("分值")
    if (length(it)) {
      if (grpl(it, "^0[7-9]-")) {
        it <- ""
      }
    } else it <- ""
  }
  if (is.null(it) || !length(it)) {
    it <- odk("score")
    if (is.null(it))
      it <- ""
  }
  it
}

od_get_info <- function(...) {
  it <- od_get(..., key = "info")
}

od_get_date <- function(file = "./mailparsed/date.md") {
  if (!file.exists(file)) {
    return(Sys.Date())
  }
  line <- readLines(file, n = 1)
  as.Date(line, "%A, %d %B %Y")
}

odb <- function(..., collapse = "") {
  n <- 0L
  res <- lapply(list(...),
    function(key) {
      res <- odk(key, if (!n) T else F)
      res
    })
  paste0(unlist(res), collapse = collapse)
}

odk <- function(key, fresh = F, file = "./mailparsed/part_1.md")
{
  if (file.exists(file)) {
    kds <- getOption("od_keywords")
    if (is.null(kds) | fresh) {
      lines <- readLines(file)
      pattern <- paste0("[a-z]+\\[\\[.*?\\]\\]")
      res <- stringr::str_extract_all(paste0(lines, collapse = " "), pattern)
      res <- unlist(res)
      res <- lapply(res,
        function(x) {
          name <- gs(x, "^([a-z]+).*", "\\1")
          value <- gs(x, ".*\\[\\[(.*?)\\]\\].*", "\\1")
          list(name, value)
        })
      kds <- lapply(res, function(x) x[[2]])
      names(kds) <- lapply(res, function(x) x[[1]])
      options(od_keywords = kds)
    }
    kds[[ key ]]
  } else {
    NULL
  }
}

revise_items <- function(data, records = NULL) {
  owd <- getwd()
  fun <- function(dirs) {
    lapply(dirs,
      function(dir) {
        setwd(dir)
        browseURL("output.pdf")
        item <- readRDS(".items.rds")
        item$lock <- T
        if (is.null(item$tags)) {
          item$tags <- ""
        }
        item <- edit(item)
        saveRDS(item, ".items.rds")
      })
  }
  tryCatch(fun(data$.dir), finally = setwd(owd))
}

# get_tags <- function(data) {
#   tags <- unlist(lapply(data$.dir,
#       function(dir) {
#         obj <- try(readRDS(paste0(dir, "/.items.rds")), T)
#         if (!inherits(obj, "try-error")) {
#           gs(strsplit(obj$tags, ",")[[1]], "\\s", "")
#         } else NULL
#       }))
#   unique(tags)
# }

od_get <- function(file = "./mailparsed/part_1.md", key = "id",
  pattern = paste0("(?<=", key, "\\{\\{).*?(?=\\}\\})"))
{
  if (file.exists(file)) {
    lines <- readLines(file)
    if (is.null(pattern)) {
      lines
    } else {
      res <- stringr::str_extract(paste0(lines, collapse = " "), pattern)
      if (is.na(res)) {
        if (key == "info") {
          return(paste0(.guess_info(lines), collapse = " "))
        } else if (key == "id") {
          return(.guess_id(lines))
        } else {
          return("")
        }
      } else {
        return(res)
      }
    }
  } else {
    return()
  }
}

.guess_client <- function(lines, excludes = c("生信", "新单", "标书"))
{
  its <- lines[ grpl(lines, "主题|Subject") ]
  its <- gs(its, ".*Subject|.*主题", "")
  maybes <- unlist(stringr::str_extract_all(its, "[\u4e00-\u9fa5]+"))
  maybes <- unique(maybes)
  maybes <- maybes[ nchar(maybes) >= 2 & nchar(maybes) <= 4 ]
  maybes <- maybes[ !grpl(maybes, paste0(excludes, collapse = "|")) ]
  if (length(maybes) > 1) {
    which <- menu(c("None", maybes), title = "Which is the client name?")
    if (which == 1) {
      client <- NULL
    } else {
      client <- maybes[ which - 1 ]
    }
  } else {
    client <- maybes
  }
  client
}

.guess_id <- function(lines) {
  res <- stringr::str_extract_all(lines, "[A-Z]{0,2}[0-9]{6,10}")
  res <- unlist(res)
  possibles <- res <- res[ !is.na(res) ]
  if (length(res)) {
    stat <- table(res)
    id <- names(stat[ stat == max(stat) ])
    if (length(id) > 1) {
      id <- id[1]
    }
    if (strx(id, "[0-9]{2}") != "20") {
      isThat <- usethis::ui_yeah(paste0("Detected ID: ", id, "\n\tUsethis?"))
      if (!isThat) {
        possibles <- c(possibles, "")
        which <- menu(possibles, title = "Use which?")
        possibles[ which ]
      } else id
    } else {
      id
    }
  } else {
    ""
  }
}

.guess_info <- function(lines) {
  raw_numlst <- as.integer(strx(lines, "^[0-9]+"))
  numlst <- raw_numlst[ !is.na(raw_numlst) ]
  if (length(numlst)) {
    isThat <- numlst[ -length(numlst) ] + 1 == numlst[-1]
    if (!any(isThat)) {
      return("")
    }
    if (table(isThat)[[ "TRUE" ]] > 3) {
      n <- 0L
      j <- 0L
      for (i in isThat) {
        j <- j + 1L
        if (n > 3) {
          break
        } else if (i) {
          n <- n + 1L
        } else {
          n <- 0L
        }
      }
      start <- which(!is.na(raw_numlst))[ j - 4 ]
      numlst <- ifelse(is.na(raw_numlst), 0, raw_numlst)
      maybeEnd <- which(numlst == max(numlst))[1]
      n <- maybeEnd + 1L
      cont <- lines[ n ]
      while (cont != "") {
        n <- n + 1L
        cont <- lines[ n ]
      }
      end <- n
      return(lines[ start:end ])
    } else {
      return("")
    }
  } else {
    return("")
  }
}

deparse_mail <- function(dir = "mail",
  savedir = "mailparsed", attsdir = "order_material",
  force = F)
{
  if (!dir.exists(dir)) {
    message("No mail directory found.")
    return()
  }
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  if (!force & file.exists(sig <- paste0(savedir, "/", ".parsed.md"))) {
    return()
  }
  testfile <- list.files(dir, full.names = T, recursive = T)[1]
  fewLines <- readLines(testfile, n = 100)
  if (!any(grpl(fewLines, "Content-Transfer-Encoding"))) {
    message("Maybe this mail is not the raw file.")
    return()
  }
  ## import python package
  bt <- e(reticulate::import_builtins())
  m <- e(reticulate::import("mailbox"))
  ## load mail
  mdir <- m$Maildir(dir)
  obj <- mdir$get_message(mdir$keys()[1]) 
  ## payload
  contents <- obj$get_payload()
  isMulti <- vapply(contents, function(x) x$is_multipart(), logical(1))
  ## all attachments
  atts <- contents[ !isMulti ]
  if (!dir.exists(attsdir)) {
    dir.create(attsdir)
  }
  if (length(atts) == 0) {
    writeLines("", paste0(attsdir, "/empty.txt"))
  } else {
    n <- 0L
    lapply(atts,
      function(att) {
        bins <- att$get_payload(decode = T)
        filename <- att$get_filename()
        if (length(filename) == 0) {
          return()
        }
        if (grpl(filename, "^=\\?UTF-8")) {
          path <- paste0(savedir, "/", filename)
        } else {
          n <<- n + 1L
          suffix <- stringr::str_extract(filename, "\\.[a-zA-Z0-9]*$")
          if (is.na(suffix)) {
            suffix <- ""
          }
          path <- paste0(attsdir, "/", "file_", n, suffix)
        }
        fp <- bt$open(path, "wb")
        fp$write(bins)
        fp$close()
      })
  }
  if (!length(list.files(attsdir))) {
    writeLines("", paste0(attsdir, "/empty.txt"))
  }
  ## multipart
  main <- contents[ isMulti ]
  if (length(main) == 0) {
    # return()
    main <- contents
  } else {
    main <- main[[1]]
    main <- main$get_payload()
  }
  n <- 0L
  files_multipart <- lapply(main,
    function(part) {
      n <<- n + 1L
      if (part$get_content_type() == "multipart/alternative") {
        part <- part$get_payload()[[1]]
      }
      bin <- part$get_payload(decode = T)
      if (part$get_content_type() == "text/plain") {
        filename <- paste0("part_", n, ".md")
      } else if (part$get_content_type() == "text/html") {
        filename <- paste0("part_", n, ".html")
      } else {
        return()
      }
      fp <- bt$open(file <- paste0(savedir, "/", filename), "wb")
      fp$write(bin)
      fp$close()
      return(file)
    })
  files_multipart <- lst_clear0(files_multipart)
  writeLines("", sig)
  ## basic information
  date <- lapply(obj$items(),
    function(lst) {
      if (lst[[ 1 ]] == "Date") {
        return(lst[[ 2 ]])
      }
    })
  date <- unlist(date)
  writeLines(date, paste0(savedir, "/", "date.md"))
}


