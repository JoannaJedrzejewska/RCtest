# =============================================================================
# Table 8 supporting audit: full exported-function inventories plus keyword
# matches for packages compared in the manuscript.
# The inventories were also reviewed manually.
# =============================================================================

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

packages_to_audit <- c(
  "RCtest",
  "forecast",
  "MCS",
  "esback",
  "scoringRules",
  "ForeComp",
  "GAS",
  "ExactVaRTest"
)

keywords <- c(
  "reality",
  "white",
  "spa",
  "superior",
  "conditional",
  "cpa",
  "diebold",
  "mariano",
  "dm",
  "crps",
  "score",
  "log",
  "klic",
  "density",
  "zp",
  "var",
  "kupiec",
  "christoffersen",
  "quantile",
  "backtest",
  "expected",
  "shortfall",
  "es",
  "model",
  "confidence",
  "mcs"
)

get_exported_function_names <- function(package_name) {
  exports <- getNamespaceExports(package_name)
  
  exports[
    vapply(
      exports,
      function(object_name) {
        exported_object <- tryCatch(
          getExportedValue(package_name, object_name),
          error = function(error) NULL
        )
        
        is.function(exported_object)
      },
      logical(1)
    )
  ]
}

inventory_rows <- list()
keyword_rows <- list()
version_rows <- list()

for (package_name in packages_to_audit) {
  installed <- requireNamespace(package_name, quietly = TRUE)
  
  if (!installed) {
    version_rows[[package_name]] <- data.frame(
      Package = package_name,
      Installed = FALSE,
      Version = NA_character_,
      Audit_Date = as.character(Sys.Date()),
      R_Version = R.version.string,
      stringsAsFactors = FALSE
    )
    
    inventory_rows[[package_name]] <- data.frame(
      Package = package_name,
      Installed = FALSE,
      Version = NA_character_,
      Function = NA_character_,
      Audit_Date = as.character(Sys.Date()),
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  package_version <- as.character(utils::packageVersion(package_name))
  
  version_rows[[package_name]] <- data.frame(
    Package = package_name,
    Installed = TRUE,
    Version = package_version,
    Audit_Date = as.character(Sys.Date()),
    R_Version = R.version.string,
    stringsAsFactors = FALSE
  )
  
  exported_functions <- sort(get_exported_function_names(package_name))
  
  inventory_rows[[package_name]] <- data.frame(
    Package = package_name,
    Installed = TRUE,
    Version = package_version,
    Function = exported_functions,
    Audit_Date = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
  
  for (keyword in keywords) {
    function_matches <- exported_functions[
      grepl(keyword, exported_functions, ignore.case = TRUE)
    ]
    
    help_result <- tryCatch(
      utils::help.search(
        pattern = keyword,
        package = package_name,
        fields = c("alias", "concept", "title"),
        agrep = FALSE
      ),
      error = function(error) NULL
    )
    
    help_matches <- character(0)
    
    if (!is.null(help_result) &&
        !is.null(help_result$matches) &&
        nrow(help_result$matches) > 0L) {
      help_matches <- unique(
        paste(
          help_result$matches[, "Package"],
          help_result$matches[, "Topic"],
          sep = "::"
        )
      )
    }
    
    keyword_rows[[paste(package_name, keyword, sep = "_")]] <- data.frame(
      Package = package_name,
      Version = package_version,
      Keyword = keyword,
      Export_Matches = if (length(function_matches) == 0L) {
        ""
      } else {
        paste(function_matches, collapse = "; ")
      },
      Help_Matches = if (length(help_matches) == 0L) {
        ""
      } else {
        paste(sort(help_matches), collapse = "; ")
      },
      Audit_Date = as.character(Sys.Date()),
      stringsAsFactors = FALSE
    )
  }
}

function_inventory <- do.call(rbind, inventory_rows)
keyword_matches <- do.call(rbind, keyword_rows)
package_versions <- do.call(rbind, version_rows)

write.csv(
  function_inventory,
  file.path(output_dir, "table6_package_function_inventory.csv"),
  row.names = FALSE
)

write.csv(
  keyword_matches,
  file.path(output_dir, "table6_package_keyword_matches.csv"),
  row.names = FALSE
)

write.csv(
  package_versions,
  file.path(output_dir, "table6_package_versions.csv"),
  row.names = FALSE
)
