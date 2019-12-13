#### ImmGen Shiny Functions ####

# this script can be sourced from the shiny app and will store all funcitons needed
# for plotting


############### EXPRESSION BARPLOTS ######################
# barplot for an individual miRNA


mirnaBarplot <- function(exp, mirna, cell_types, cell_classes){
  if(is.null(cell_types) && is.null(cell_classes)){
    ggplot() + geom_blank() + theme_few()
  }
  else{
    # extract sample names that are either checked or are checked via class
    if(!is.null(cell_classes)){
    cell_table <- distinct(exp, cell_type, cell_class) %>% collect()
    cell_types_join <- left_join(data.frame(cell_class = cell_classes), cell_table,
                                 by = "cell_class") %>% distinct(cell_type) %>%
       .$cell_type
    cell_types_join <- unique(c(cell_types_join, cell_types))
    }
    else{
      cell_types_join <- cell_types
    }
  exp %>%
    filter(miRNA == mirna) %>% collect() %>%
    left_join(data.frame(cell_type = cell_types_join), ., by = "cell_type") %>%
    mutate(Sample = fct_relevel(factor(Sample), as.character(sample_order[sample_order %in% Sample]))) %>%
    ggplot(., aes(x = Sample, y = expression)) +
    geom_bar(aes(fill = cell_type), color = 'black', position = 'dodge', stat = 'identity') +
    theme_few() +
    labs(title = mirna, x = "", y = "") +
    theme(plot.title = element_text(hjust = .5, size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14)) +
    guides(fill = FALSE)

  }
}

############ TRANSCRIPT VIEW PLOTS ##################

# this will have to be edited, but right now I am using a semi join to filter
# down the number of things I'm plotting to just the family names that match

## family names need to be worked out, there seems to be a discrepancy

## new version of transcript repression barplot
### will use dplyr for now but may want to switch to data.table for speed
plotGeneTscan <- function(rep_tscript, tab, utr,
                          exp_diff_check = F,
                          exp_diff_thresh = 3,
                          context_thresh = 0,
                          exp_thresh = 0, cell_types_chosen = NULL){

  tscript_length <- dplyr::filter(utr, Transcript_ID == rep_tscript) %>% collect() %>% .$length.calc

  g <- ggplot(na.omit(tab), aes(x = UTR_start, y = mean_expression)) +
    xlim(c(0, tscript_length)) +
    ylim(c(0, max(tab$mean_expression))) +
    geom_bar(position = "identity", stat = "identity",
             width = 7, aes(fill = `context++_score` * -1)) +
    geom_text_repel(data = filter(na.omit(tab), mean_expression > exp_thresh),
                    aes(label = miRNA), size = 3, stat = "identity", ylim = c(0,NA)) +
    geom_point(shape = 21, fill = 'white') +
    labs(x = "Base Pairs from Start of 3'UTR", y = "Log2 Abundance") +
    theme_few() +
    geom_hline(yintercept = 0) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = .5),
          legend.position = "bottom",
          strip.text = element_text(size = 24),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 16)
    ) +
    guides(color = guide_legend(override.aes = list(size=6))) +
    scale_fill_gradient(low = 'grey90', high = 'grey0', name = "context++ score") +
    facet_wrap(~ cell_type, ncol = 1)

  if(exp_diff_check == T){

    # filter the targets for ones with expression differences
    tab.filt <- filter(na.omit(tab), abs(exp_diff) > exp_diff_thresh) %>%
      mutate(diff_dir = as.character(sign(exp_diff)))
    # get the correct names for the direction of expression change

    color_labels <- paste(cell_types_chosen, "up")
    names(color_labels) <- c('-1', '1')

    g <- g + geom_point(data = tab.filt,
                   aes(color = diff_dir), shape = 20) +
      scale_color_manual(name = "",
                         values = c('dodgerblue', 'red'),
                         #breaks = c('-1', '1'),
                         labels = color_labels)
  }
  g
}


##### Transcript Repression Barplot ####
## Work in progress, probably will not use
# indTranscriptBPshiny <- function(transcript, cell_class, tscript_mirna_reg){
#   tscript_mirna_reg %>%
#     filter(Transcript_ID == transcript & grepl(paste(get(cell_class), collapse = "|"), cell_type)) %>%
#     ggplot(aes(x = cell_type, y = tscript_repression)) +
#     geom_bar(stat = "identity", aes(fill = cell_type)) +
#     labs(title = paste(transcript, "repression in", cell_class, sep = " "), x = "", y = "Summed Linear Expression") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     guides(fill = FALSE)
# }
#
# indTranscriptBPshiny("ENSMUST00000028111.4", "T_cells", tscript_mirna_singSum)
