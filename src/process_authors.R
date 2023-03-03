process_authors <- function(author_file, with_degrees = FALSE) {
    authors <- yaml::read_yaml(author_file)

    affiliations <- 1:length(unique(unlist(lapply(
        authors$author, function(x) {
            x$affiliations
        }
    ))))

    names(affiliations) <- unique(unlist(lapply(
        authors$author, function(x) {
            x$affiliations
        }
    )))

    for (i in 1:length(authors$author)) {
        authors$author[[i]]$superscript <-
            paste0(
                affiliations[authors$author[[i]]$affiliations],
                collapse = ","
            )
        authors$author[[i]]$label <-
            paste0(
                authors$author[[i]]$name, "^", authors$author[[i]]$superscript, "^"
            )
        if (with_degrees) {
            authors$author[[i]]$label <-
                paste0(
                    authors$author[[i]]$label, ", ", paste0(authors$author[[i]]$degree, collapse = ", ")
                )
        }
        if ("corresponding" %in% authors$author[[i]]$attributes) {
            authors$author[[i]]$label <-
                paste0(
                    authors$author[[i]]$name, "^", authors$author[[i]]$superscript, "\\*^"
                )
        }
    }

    parsed_authors <- paste0(unlist(lapply(
        authors$author, function(x) {
            x$label
        }
    )), collapse = ", ")

    for (i in names(authors$affiliations)) {
        authors$affiliations[[i]]$number <- unname(affiliations[i])
        authors$affiliations[[i]]$label <-
            paste0(
                "**",
                authors$affiliations[[i]]$number,
                "** ",
                str_c(
                    authors$affiliations[[i]]$department,
                    authors$affiliations[[i]]$name,
                    authors$affiliations[[i]]$city,
                    authors$affiliations[[i]]$state,
                    authors$affiliations[[i]]$country,
                    sep = ", "
                )
            )
    }


    parsed_affiliations <- paste0(unlist(lapply(
        names(affiliations),
        function(x) {
            authors$affiliations[[x]]$label
        }
    )),
    collapse = "  \n"
    )

    abstract <- authors$abstract
    acknowledgements <- authors$acknowledgements
    figure_legends <- authors$figure_legends
    cois <- authors$cois
    return(list(
        parsed_authors = parsed_authors,
        parsed_affiliations = parsed_affiliations,
        abstract = abstract,
        acknowledgements = acknowledgements,
        figure_legends = figure_legends,
        cois = cois
    ))
}
