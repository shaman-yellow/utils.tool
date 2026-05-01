# ==========================================================================
# workflow of dynamic
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_dynamic <- setClass("job_dynamic", 
  contains = c("job"),
  prototype = prototype(
    pg = "dynamic",
    info = c("https://www.gromacs.org/"),
    cite = "",
    method = "",
    tag = "dynamic",
    analysis = "分子动力学模拟"
    ))

setGeneric("asjob_dynamic",
  function(x, ...) standardGeneric("asjob_dynamic"))

setMethod("asjob_dynamic", signature = c(x = "list"),
  function(x){
    if (!is.null(names(x))) {
      x <- c(unname(x), sigs = list(names(x)))
    }
    do.call(asjob_dynamic, x)
  })

setMethod("asjob_dynamic", signature = c(x = "job_vina"),
  function(x, ..., sigs = NULL)
  {
    if (...length()) {
      objects <- list(x, ...)
    } else {
      objects <- list(x)
    }
    meta <- lapply(objects,
      function(obj) {
        if (!is(obj, "job_vina")) {
          stop('!is(obj, "job_vina").')
        }
        if (obj@step < 8L || is.null(obj$res_dock_merge)) {
          stop('obj@step < 8L || is.null(obj$res_dock_merge).')
        }
        obj$res_dock_merge
      })
    meta <- dplyr::bind_rows(meta, .id = ".id")
    if (is.null(sigs)) {
      sigs <- vapply(
        objects, FUN.VALUE = character(1),
        function(x) {
          sig(x)
        }
      )
    }
    meta$sig <- rep(sigs, lengths(split(meta$.id, meta$.id)))
    x <- .job_dynamic()
    x$metadata <- tibble::as_tibble(meta)
    x <- methodAdd(x, "为了进一步了解化合物与蛋白之间的作用机制，以分子动力学模拟来评估化合物与蛋白复合物的稳定性以及动力学特性。将 AutoDock vina 的分子对接结果文件格式化为 GROMACS 软件输入。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_dynamic"),
  function(x){
    step_message("Prepare your data with function `asjob_dynamic`.")
  })

setMethod("step1", signature = c(x = "job_dynamic"),
  function(x) {
    step_message("Do nothing.")
    x <- methodAdd(x, "使用 GROMACS 2024.4 软件进行分子动力学模拟。体系采用 AMBER99SB-ILDN 力场，并使用 TIP3P 水模型构建溶剂环境。构建溶剂盒时确保蛋白质与盒边之间的最小距离为 1 nm，并通过添加适量离子以维持体系电中性。首先采用最速下降法对体系进行能量最小化。随后依次进行等温等体积（NVT）和等温等压（NPT）平衡过程，各阶段约持续 100 ps。在平衡过程中，温度耦合采用 V-rescale 方法，参考温度设定为 300 K；压力耦合采用 Parrinello–Rahman 方法。积分步长设定为 2 fs，并采用 LINCS 算法约束键长。在体系达到平衡后，进行 20 ns 的分子动力学模拟。基于模拟轨迹计算蛋白–配体复合物的 RMSD、蛋白 RMSF、体系总能量以及氢键数量等参数。")
    return(x)
  })


setMethod("step2", signature = c(x = "job_dynamic"),
  function(x, dir_results, fun_detect = .detect_gromacs_xvg_project.hb, ...)
  {
    step_message("Load results and plot.")
    x$metadata <- fun_detect(x$metadata, dir_results)
    x$metadata <- dplyr::mutate(
      x$metadata,
      short_name = stringr::str_trunc(Ingredient_name, 30),
      short_pdb_id = ifelse(
        grpl(PDB_ID, "^AF-", TRUE), strx(PDB_ID, "[Aa][Ff]-[^-]+"), PDB_ID
      ),
      labels = glue::glue("{hgnc_symbol} ({short_pdb_id}) + {short_name} (CID: {PubChem_id})"),
      labels = stringr::str_wrap(labels, 30)
    )
    x$md_data <- .load_as_md_data(
      x$metadata$dir_dynamic, x$metadata$labels
    )
    ps <- .plot_alls_md_data(x$md_data, ...)
    p.rmsd <- set_lab_legend(
      ps$p.rmsd,
      glue::glue("{x@sig} gromacs RMSD plot"),
      glue::glue(
        "蛋白–配体复合物RMSD变化曲线|||展示各复合物在分子动力学模拟过程中整体构象偏移随时间的变化情况。横坐标为模拟时间（ns），纵坐标为RMSD（nm），不同颜色代表不同复合物体系。该图用于评估体系在模拟过程中的构象稳定性及收敛情况，RMSD波动较小且趋于平台通常表示体系较为稳定，而持续波动或上升则提示可能存在构象调整或稳定性较差。"
      )
    )
    p.rmsf <- set_lab_legend(
      ps$p.rmsf,
      glue::glue("{x@sig} gromacs RMSF plot"),
      glue::glue(
        "蛋白残基RMSF波动分布图|||展示各复合物体系中蛋白质各残基在模拟过程中的柔性分布情况。横坐标为残基编号（按各体系重新编号），纵坐标为RMSF（nm），不同颜色代表不同复合物。该图用于分析蛋白不同区域的动态特性，较高的RMSF值通常对应柔性较大的区域（如环区或末端），而较低的RMSF则表明结构相对刚性稳定。"
      )
    )
    p.energy <- set_lab_legend(
      ps$p.energy,
      glue::glue("{x@sig} gromacs energy plot"),
      glue::glue(
        "体系总能量变化曲线|||展示各复合物体系在分子动力学模拟过程中总能量随时间的变化趋势。横坐标为模拟时间（ns），纵坐标为能量（原始能量经组内中心化后的ΔEnergy）。图中每个子图对应一个复合物体系，并在左下角标注其平均 (原始) 能量及标准差。该图用于评估体系的热力学稳定性及能量收敛情况，能量波动较小通常表示体系趋于稳定，而较大波动可能提示尚未充分平衡或存在构象变化。"
      )
    )
    p.hbond <- set_lab_legend(
      ps$p.hbond,
      glue::glue("{x@sig} gromacs hydrogen bond plot"),
      glue::glue(
        "氢键数量变化曲线|||展示各蛋白–配体复合物在模拟过程中氢键数量随时间的变化情况。横坐标为模拟时间（ns），纵坐标为氢键数量，不同子图对应不同复合物体系。该图用于评估配体与蛋白之间相互作用的稳定性，氢键数量较高且波动较小通常表明结合较为稳定，而氢键数量较低或波动较大则提示相互作用较弱或不稳定。"
      )
    )
    snap <- .stat_gromacs_md_data(
      x$md_data, group = x$metadata$sig,
      p.rmsd = p.rmsd,
      p.rmsf = p.rmsf,
      p.hbond = p.hbond,
      p.energy = p.energy
    )
    x <- plotsAdd(x, p.rmsd, p.rmsf, p.hbond, p.energy)
    x <- snapAdd(x, "{snap}")
    return(x)
  })

.detect_gromacs_xvg_project.hb <- function(metadata, 
  dir, pattern = "/analysis_", col_match = "Combn")
{
  dirs <- list.dirs(dir, full.names = TRUE, recursive = TRUE)
  dirs <- grpf(dirs, pattern)
  if (length(dirs) != nrow(metadata)) {
    stop('length(dirs) != nrow(metadata), results number not matched.')
  }
  lapply(dirs,
    function(dir) {
      files <- list.files(dir, "xvg$")
      if (length(files) < 4L) {
        stop(glue::glue("length(files) < 4L, directory is {dir}"))
      }
    })
  message(glue::glue("Match pairs ID with directory name:"))
  matches <- match_strings(
    metadata[[ col_match ]], basename(dirs),
    max_dist = .5, method = "cosine", keep_identical = TRUE
  )
  if (any(is.na(matches$matched_index_in_y))) {
    stop('any(is.na(matches$matched_index_in_y)).')
  }
  metadata$dir_dynamic <- dirs[ matches$matched_index_in_y ]
  metadata
}

.plot_alls_md_data <- function(obj_md, lst_color = NULL, ...)
{
  if (is.null(lst_color)) {
    lst_color <- color_set()[seq_along(obj_md$group)]
  }
  layout <- wrap_layout(NULL, length(obj_md$group))
  lst <- list(
    p.rmsd = .plot_rmsd_md_data(obj_md, lst_color),
    p.rmsf = .plot_rmsf_md_data(obj_md, lst_color),
    p.hbond = .plot_hbond_md_data(obj_md, lst_color, ncol = layout$ncol),
    p.energy = .plot_energy_md_data(obj_md, lst_color, ncol = layout$ncol)
  )
  lst[1:2] <- lapply(lst[1:2], wrap)
  lst[3:4] <- lapply(lst[3:4],
    function(x) {
      add(layout, x)
    })
  return(lst)
}

# ==============================
# S3 loader
# ==============================

.load_as_md_data <- function(
  vec_path,
  vec_group,
  file_rmsd   = "rmsd.xvg",
  file_rmsf   = "rmsf.xvg",
  file_energy = "energy_total.xvg",
  file_hbond  = "hbond_num.xvg"
)
{

  if (length(vec_path) != length(vec_group)) {
    stop("vec_path and vec_group must have same length")
  }

  message(glue::glue("Start loading MD data: {length(vec_path)} groups"))

  .load_one <- function(file_name, fun_read) {

    lst_data <- lapply(
      seq_along(vec_path),
      function(i) {

        path_now <- file.path(vec_path[i], file_name)

        if (!file.exists(path_now)) {
          message(glue::glue("Missing file: {path_now}"))
          return(NULL)
        }

        data_now <- fun_read(path_now)
        data_now$group <- vec_group[i]

        return(data_now)
      }
    )

    lst_data <- Filter(Negate(is.null), lst_data)

    if (length(lst_data) == 0L) {
      message(glue::glue("No valid file for {file_name}"))
      return(NULL)
    }

    data_all <- do.call(rbind, lst_data)

    message(glue::glue("{file_name}: {nrow(data_all)} rows"))

    return(data_all)
  }

  lst_all <- list(
    rmsd   = .load_one(file_rmsd,   .read_xvg),
    rmsf   = .load_one(file_rmsf,   .read_xvg),
    energy = .load_one(file_energy, .read_xvg),
    hbond  = .load_one(file_hbond,  .read_xvg)
  )

  structure(
    list(
      data  = lst_all,
      group = vec_group,
      path  = vec_path
    ),
    class = "md_data"
  )
}


# ==============================
# checker
# ==============================

.check_md_data <- function(obj_md) {

  if (!is(obj_md, "md_data")) {
    stop("Input must be md_data object")
  }

  return(TRUE)
}


# ==============================
# RMSD
# ==============================

.plot_rmsd_md_data <- function(
  obj_md,
  lst_color, ...
) {

  .check_md_data(obj_md)

  data_all <- obj_md$data$rmsd

  if (is.null(data_all)) {
    stop("No RMSD data available")
  }

  p <- ggplot(
    data_all,
    aes(x = `Time (ns)`, y = `RMSD (nm)`, color = group)
  ) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    # facet_wrap(~ group, ...) +
    scale_color_manual(values = lst_color) +
    labs(x = "Time (ns)", y = "RMSD (nm)", color = "Group") +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)
}


# ==============================
# RMSF
# ==============================

.plot_rmsf_md_data <- function(
  obj_md,
  lst_color, ...
) {

  .check_md_data(obj_md)

  data_all <- obj_md$data$rmsf

  if (is.null(data_all)) {
    stop("No RMSF data available")
  }

  data_all$Residue <- ave(
    seq_len(nrow(data_all)),
    data_all$group,
    FUN = seq_along
  )

  p <- ggplot(
    data_all,
    aes(x = Residue, y = `(nm)`, color = group)
  ) +
    geom_line(linewidth = 0.7) +
    # facet_wrap(~ group, ...) +
    scale_color_manual(values = lst_color) +
    labs(x = "Residue", y = "RMSF (nm)", color = "Group") +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)
}

.plot_energy_md_data <- function(
  obj_md,
  lst_color,
  use_center = TRUE,
  use_smooth = TRUE,
  n_roll = 100L, ...
) {

  .check_md_data(obj_md)

  data_all <- obj_md$data$energy

  if (is.null(data_all)) {
    stop("No energy data available")
  }

  if ("Time (ps)" %in% colnames(data_all)) {
    data_all$`Time (ps)` <- data_all$`Time (ps)` / 1000
    colnames(data_all)[ colnames(data_all) == "Time (ps)" ] <- "Time (ns)"
  }

  if (!"Total Energy" %in% colnames(data_all)) {
    stop("Column 'Total Energy' not found")
  }

  message(glue::glue(
    "Energy plot | center = {use_center}, smooth = {use_smooth}, n_roll = {n_roll}"
  ))

  vec_energy <- data_all$`Total Energy`
  vec_group  <- data_all$group

  # ---------- compute mean + sd (use ORIGINAL energy) ----------
  data_stat <- stats::aggregate(
    data_all$`Total Energy`,
    by = list(group = data_all$group),
    FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
  )

  data_stat <- data.frame(
    group = data_stat$group,
    mean  = data_stat$x[, "mean"],
    sd    = data_stat$x[, "sd"]
  )

  # ---------- center ----------
  if (use_center) {
    vec_energy <- vec_energy -
      ave(vec_energy, vec_group, FUN = mean)
  }

  # ---------- smooth ----------
  if (use_smooth) {

    if (!requireNamespace("zoo", quietly = TRUE)) {
      message("Package 'zoo' not found, skip smoothing")
    } else {

      vec_energy <- ave(
        vec_energy,
        vec_group,
        FUN = function(x) {

          if (length(x) < n_roll) {
            return(x)
          }

          zoo::rollmean(
            x,
            k = n_roll,
            fill = "extend"
          )
        }
      )
    }
  }

  data_all$energy_plot <- vec_energy

  # ---------- label position (per facet) ----------
  data_label <- stats::aggregate(
    data_all$energy_plot,
    by = list(group = data_all$group),
    FUN = function(x) {
      x[ round(length(x) * 0.9) ]
    }
  )

  colnames(data_label)[2] <- "y"

  # x position
  data_label$x <- stats::aggregate(
    data_all$`Time (ns)`,
    by = list(group = data_all$group),
    FUN = function(x) max(x) * 0.95
  )$x

  # merge stats
  data_label <- merge(data_label, data_stat, by = "group")

  # label text
  data_label$label <- sprintf(
    "mean = %.2f\nsd = %.2f",
    data_label$mean,
    data_label$sd
  )

  # ---------- y label ----------
  str_ylabel <- "Energy (kcal/mol)"

  if (use_center && use_smooth) {
    str_ylabel <- glue::glue(
      "ΔEnergy (kcal/mol)\n(mean-centered within group, {n_roll}-point rolling mean)"
    )
  } else if (use_center) {
    str_ylabel <- "ΔEnergy (kcal/mol)\n(mean-centered within group)"
  } else if (use_smooth) {
    str_ylabel <- glue::glue(
      "Energy (kcal/mol)\n({n_roll}-point rolling mean)"
    )
  }

  p <- ggplot(
    data_all,
    aes(
      x = `Time (ns)`,
      y = energy_plot,
      color = group
    )
  ) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_text(
      data = data_label,
      aes(
        x = -Inf,
        y = -Inf,
        label = label
        ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = -0.5,
      size = 3
    ) +
    scale_color_manual(values = lst_color) +
    labs(
      x = "Time (ns)",
      y = str_ylabel,
      color = "Group"
    ) +
    facet_wrap(~ group, ...) +
    theme_minimal() +
    theme(legend.position = "none")

  return(p)
}

# ==============================
# HBond
# ==============================

.plot_hbond_md_data <- function(
  obj_md,
  lst_color, ...
) {

  .check_md_data(obj_md)

  data_all <- obj_md$data$hbond

  if (is.null(data_all)) {
    stop("No hbond data available")
  }

  if ("Time (ps)" %in% colnames(data_all)) {
    data_all$`Time (ps)` <- data_all$`Time (ps)` / 1000
    colnames(data_all)[ colnames(data_all) == "Time (ps)" ] <- "Time (ns)"
  }

  p <- ggplot(
    data_all,
    aes(x = `Time (ns)`, y = `Hydrogen bonds`, color = group)
  ) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ group, ...) +
    scale_color_manual(values = lst_color) +
    labs(x = "Time (ns)", y = "HBonds", color = "Group") +
    theme_minimal() +
    theme(legend.position = "none")

  return(p)
}

.stat_gromacs_md_data <- function(
  obj_md,
  labels = obj_md$group,
  group = NULL,
  p.rmsd = NULL,
  p.rmsf = NULL,
  p.hbond = NULL,
  p.energy = NULL,
  bot = FALSE,
  sum = FALSE
)
{

  .check_md_data(obj_md)

  lst_data <- obj_md$data

  data_rmsd   <- lst_data$rmsd
  data_rmsf   <- lst_data$rmsf
  data_hbond  <- lst_data$hbond
  data_energy <- lst_data$energy

  if (is.null(group)) {
    group <- rep("All", length(unique(data_rmsd$group)))
  }

  vec_group <- unique(data_rmsd$group)
  data_map <- data.frame(
    group = vec_group,
    label = labels,
    super = group
  )

  # ============================
  # column detection
  # ============================

  col_rmsd  <- colnames(data_rmsd)[2L]
  col_rmsf  <- colnames(data_rmsf)[2L]
  col_hbond <- colnames(data_hbond)[2L]

  # ============================
  # summarise
  # ============================

  data_rmsd_stat <- dplyr::summarise(
    dplyr::group_by(data_rmsd, group),
    rmsd_mean = mean(.data[[col_rmsd]]),
    rmsd_sd   = stats::sd(.data[[col_rmsd]]),
    .groups = "drop"
  )

  data_rmsf_stat <- dplyr::summarise(
    dplyr::group_by(data_rmsf, group),
    rmsf_mean = mean(.data[[col_rmsf]]),
    .groups = "drop"
  )

  data_hbond_stat <- dplyr::summarise(
    dplyr::group_by(data_hbond, group),
    hbond_mean = mean(.data[[col_hbond]]),
    .groups = "drop"
  )

  data_energy_stat <- dplyr::summarise(
    dplyr::group_by(data_energy, group),
    energy_sd = stats::sd(`Total Energy`),
    .groups = "drop"
  )

  data_stat <- Reduce(
    function(x, y) merge(x, y, by = "group"),
    list(data_rmsd_stat, data_rmsf_stat, data_hbond_stat, data_energy_stat)
  )

  data_stat <- merge(data_stat, data_map, by = "group")

  # ============================
  # scoring
  # ============================

  rank_scale <- function(x) rank(x) / length(x)

  data_stat$score <-
    (1 - rank_scale(data_stat$rmsd_sd)) * 0.35 +
    (1 - rank_scale(data_stat$rmsf_mean)) * 0.2 +
    rank_scale(data_stat$hbond_mean)     * 0.3 +
    (1 - rank_scale(data_stat$energy_sd))* 0.15

  # ============================
  # reference
  # ============================

  ref_all <- aref(list(p.rmsd, p.rmsf, p.hbond, p.energy))

  # ============================
  # main description（核心改进）
  # ============================

  txt_main <- glue::glue(
    "基于分子动力学模拟结果{ref_all}，对各蛋白–配体复合物的构象稳定性与相互作用特征进行综合评估。\
其中，RMSD用于衡量整体构象偏移程度及收敛性，较低的RMSD波动（标准差较小）通常表示体系构象更稳定；\
RMSF反映蛋白各残基的柔性分布，较低的平均RMSF提示结构区域偏向刚性；\
氢键数量用于表征配体与蛋白之间的相互作用强度，数量较多且稳定通常意味着结合更加紧密；\
总能量的标准差则用于评估体系的热力学波动程度，较低波动通常提示体系趋于稳定状态。"
  )

  # ============================
  # group-wise evaluation
  # ============================

  lst_super <- split(data_stat, data_stat$super)

  txt_group <- sapply(names(lst_super), function(name) {

    df <- lst_super[[ name ]]
    df <- df[order(-df$score), ]

    top <- df[1L, ]
    bot_row <- df[nrow(df), ]

    txt_top <- glue::glue(
      "{top$label} 在 {name} 组中表现较优，其RMSD波动为 {round(top$rmsd_sd,3)} nm，\
表明整体构象偏移较小；同时其平均RMSF为 {round(top$rmsf_mean,3)} nm，\
提示蛋白结构整体较为刚性；此外，其氢键平均数量为 {round(top$hbond_mean,1)}，\
反映出较强且持续的分子间相互作用；能量波动（SD）为 {round(top$energy_sd,1)}，\
说明体系在模拟过程中具有较好的热力学稳定性。"
    )

    txt_bot <- ""

    if (bot) {

      if (bot_row$score < quantile(df$score, 0.25)) {

        txt_bot <- glue::glue(
          "{bot_row$label} 在 {name} 组中相对波动略高，其RMSD标准差为 {round(bot_row$rmsd_sd,3)} nm，\
且氢键数量相对较少（{round(bot_row$hbond_mean,1)}），\
提示该体系在构象或相互作用维持方面存在一定动态变化。"
        )
      }
    }

    paste(txt_top, txt_bot, sep = " ")
  })

  txt_group <- paste(txt_group, collapse = "\n\n")

  # ============================
  # global comparison
  # ============================

  best_global <- data_stat[which.max(data_stat$score), ]

  txt_global <- glue::glue(
    "在所有体系中，{best_global$label} 综合表现最为突出，\
其在RMSD波动控制、结构刚性维持及氢键相互作用数量方面均表现较优，\
同时伴随较低的能量波动，提示其复合物在动力学尺度上具有更高的稳定性。"
  )

  # ============================
  # summary
  # ============================

  txt_summary <- ""

  if (sum) {
    txt_summary <- "综合来看，不同复合物在构象稳定性与相互作用维持能力上存在一定差异，优选体系在多个具体指标上均表现出一致优势。"
  }

  txt_all <- paste(
    txt_main,
    txt_group,
    txt_global,
    txt_summary,
    sep = "\n\n"
  )

  return(txt_all)
}

# ==============================
# XVG readers
# ==============================

.read_xvg <- function(path_file) {

  vec_lines <- readLines(path_file)

  vec_data <- vec_lines[ !grepl("^[#@]", vec_lines) ]

  # x label
  vec_x <- vec_lines[ grepl("@\\s*xaxis\\s+label", vec_lines) ]
  if (length(vec_x) > 0L) {
    name_x <- gsub('.*label\\s*"(.*?)".*', "\\1", vec_x[1L])
  } else {
    name_x <- "X"
  }

  # y label (fallback)
  vec_y <- vec_lines[ grepl("@\\s*yaxis\\s+label", vec_lines) ]
  if (length(vec_y) > 0L) {
    name_y <- gsub('.*label\\s*"(.*?)".*', "\\1", vec_y[1L])
  } else {
    name_y <- "Y"
  }

  # legends (multi-line support)
  vec_leg <- vec_lines[ grepl("@\\s*s\\d+\\s+legend", vec_lines) ]

  if (length(vec_leg) > 0L) {
    vec_y_names <- gsub('.*legend\\s*"(.*?)".*', "\\1", vec_leg)
  } else {
    vec_y_names <- name_y
  }

  data_out <- read.table(text = vec_data, header = FALSE)

  n_col <- ncol(data_out)

  # construct column names
  if (n_col == 1L) {
    vec_col <- name_x
  } else if (n_col == 2L) {
    vec_col <- c(name_x, vec_y_names[1L])
  } else {
    # multi Y columns
    if (length(vec_y_names) == (n_col - 1L)) {
      vec_col <- c(name_x, vec_y_names)
    } else {
      vec_col <- c(name_x, paste0("Y", seq_len(n_col - 1L)))
    }
  }

  colnames(data_out) <- vec_col

  return(data_out)
}

.make_color_mapping <- function(
  vec_group,
  seed = 12345L
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_group <- length(vec_group)

  if (n_group == 0L) {
    stop("vec_group is empty")
  }

  # generate evenly spaced hues with random offset
  val_offset <- stats::runif(1L, 0, 360)

  vec_hue <- (
    seq(0, 360 - 360 / n_group, length.out = n_group) + val_offset
  ) %% 360

  # shuffle for better separation in facet
  vec_hue <- sample(vec_hue)

  vec_color <- grDevices::hcl(
    h = vec_hue,
    c = 50,
    l = 80
  )

  lst_color <- stats::setNames(vec_color, vec_group)

  message(glue::glue("Generated {n_group} colors"))

  return(lst_color)
}

