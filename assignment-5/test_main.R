#!/usr/bin/Rscript
source("main.R")
library(testthat)

test_data <- as_tibble(read.csv("data/test.tsv", sep=""))
td2 <- test_data[1:1000,]

describe("read_data()", {
  data <- read_data("data/verse_counts.tsv")
  type <- unlist(data %>% summarize_all(class))
  
  it("returns a tibble, not a dataframe", {
    expect_true(is_tibble(data))
  })
  it("has no rownames, tibbles do not use rownames", {
    expect_false(has_rownames(data))
  })
  it("returns a tibble with 55416 rows and 9 columns", {
    expect_equal(dim(data), c(55416, 9))
  })
  it("has the correct column names", {
    expect_setequal(colnames(data), c("gene", "vP0_1", "vP0_2", "vP4_1", "vP4_2", "vP7_1", "vP7_2", "vAd_1", "vAd_2"))
  })
  it("has a gene column that stores characters", {
    expect_equal(type[['gene']], "character")
  })
  it("returns sample columns that are all numeric", {
    expect_true(unique(type[!names(type) %in% c('gene')]) == c('numeric'))
  })
})

describe("filter_zero_var_genes()", {
  col <- length(colnames(test_data))
  row <- 999
  
  zero_var_filtered <- filter_zero_var_genes(test_data)
  user_removed <- setdiff(test_data$gene, zero_var_filtered$gene)
  expected_removed <- c("CJNN53716", "WHGN65413")
  
  it("should return a tibble", {
    expect_true(is_tibble(zero_var_filtered))
  })
  it("should have exactly 999 rows and 5 columns", {
    expect_equal(dim(zero_var_filtered), c(999, 5))
  })
  it("if this test fails, you are filtering too many genes", {
    expect_false(length(user_removed) > 2)
  })
  it("if this test fails, you are not filtering enough genes", {
    expect_false(length(user_removed) < 2)
  })
  it("should remove both CJNN53716 and WHGN65413", {
    expect_setequal(user_removed, expected_removed)
  })
})

describe("timepoint_from_sample()", {
  it("works on the example case of vAd_1 and returns Ad", {
    expect_equal(timepoint_from_sample("vAd_1"), "Ad")
  })
  it("works on a second example case of vP7_2 and returns P7", {
    expect_equal(timepoint_from_sample("vP7_2"), "P7")
  })
  it("works on strings with numbers not ending in 1 or 2", {
    expect_equal(timepoint_from_sample("vP7_8"), "P7")
  })
  it("works on capital letters and numbers that aren't found in the samples", {
    expect_equal(timepoint_from_sample('vBd_7'), "Bd")
  })
  it("works when there is also a 1 somewhere in the string before the delimiter", {
    expect_equal(timepoint_from_sample("vJ1_1"), "J1")
  })
})
  
describe("sample_replicate()", {
  it("returns the right string in the most basic case", {
    expect_equal(sample_replicate("vAd_1"), "1")
  })
  it("correctly ignores numbers that are before the delimiter", {
    expect_equal(sample_replicate("vP4_1"), "1")
  })
  it("correctly captures any number, not just 1 or 2", {
    expect_equal(sample_replicate("vP3_2"), "2")
    expect_equal(sample_replicate("vPk_8"), "8")
  })
  # #Cases that don't have to work, correct functions may return differently depending on method
  # sample_replicate("vPk_82") # ("82", "8", "2", NULL, ...) # multiple numbers after the `_`
  # sample_replicate("v8k_8g") # ("8g", "8", "g", NA, ...) # It's okay if there are letters after `_`
  # sample_replicate("vLj_1_19") # ("1", "19", "9", ... ) # Returning "1_19" instead could be valid, just use `expect_equal(sample_replicate("vLj_1_19"), "1_19")`
  # sample_replicate("v1_An") # ("An", "A", NULL, ...) # doesn't have to be numbers after the `_` at all
  # sample_replicate("v_<3") #("<3", "3", NULL, ...) #non alphanumerical for kicks
})

describe("meta_info_from_labels()", {
  desired_output <-  tibble(sample=c("vK3_2", "vMr_1","vK3_1", "vMr_2"), 
                            timepoint=c("K3", "Mr", "K3", "Mr"), 
                            replicate=c("2","1","1","2"))
  
  function_labels <- meta_info_from_labels(c("vK3_2", "vMr_1","vK3_1", "vMr_2"))
  
  function_timepoints <- setNames(c(function_labels$timepoint), c(function_labels$sample))
  function_replicates <- setNames(c(function_labels$replicate), c(function_labels$sample))
  
  timepoint_vct <- setNames(c(desired_output$timepoint), c(desired_output$sample))
  replicate_vct <- setNames(c(desired_output$replicate), c(desired_output$sample))
  
  it("returns a tibble, not a dataframe", {
    expect_true(is_tibble(function_labels))
  })
  it("returns exactly 3 columns and 4 rows based on the data", {
    expect_equal(dim(function_labels), c(4, 3))
  })
  it("returns the right column names sample, timepoint and replicate", {
    expect_setequal(colnames(function_labels), colnames(desired_output))
  })
  it("returns the right row information parsed from the correspondine sample name", {
    expect_mapequal(timepoint_vct, function_timepoints)
    expect_mapequal(replicate_vct, function_replicates)
  })
})

describe("get_library_size()", {
  expected1 <- tibble("vK3_1" = 254130, 
                      "vMr_2" = 284628,
                      "vK3_2" = 267519, 
                      "vMr_1" = 283515 )
  
  expected2 <- tibble(sample = c("vK3_1", "vMr_2", "vK3_2", "vMr_1"),
                      value = c(254130,284628,267519,283515))
  
  lib_size <- get_library_size(test_data)
  dim_lib_size <- dim(lib_size)
  
  it("returns a tibble, not a dataframe and not a named vector", {
    expect_true(is_tibble(lib_size))
  })
  
  it("returns tibbles with either 4 rows and 2 columns, or 1 row and 4 columns", {
    expect_true(xor(isTRUE(all.equal(dim_lib_size, c(1, 4))), isTRUE(all.equal(dim_lib_size, c(4, 2)))))
  })
  
  if(isTRUE(all.equal(dim_lib_size, c(1, 4)))) {
    expected_vct <- expected1 %>% dplyr::slice(1) %>% unlist()
    function_vct <- lib_size %>% dplyr::slice(1) %>% unlist()
    
    it("returns the correct summed values for each column if tibble 1x4", {
      expect_mapequal(expected_vct, function_vct)
    })
    
  } else if(isTRUE(all.equal(dim_lib_size, c(4, 2)))) {
    expected_vct <- expected2 %>% pull(value, sample)
    function_vct <- lib_size %>% pull(value, sample)
    
    it("returns the correct summed values for each column if tibble 4x2", {
      expect_mapequal(expected_vct, function_vct)
    })
  }
})

describe("normalize_by_cpm()", {
  cpm <- normalize_by_cpm(td2)%>% 
    mutate(across(where(is.numeric), \(x) round(x, 3)))
  
  sample_answers <-tibble(gene = c("SZUG00090", "JIZU74666", "AWDL96666", "HDLU58875"), 
                          vK3_2 = c(44.8566270059323, 59.8088360079097, 59.8088360079097, 112.141567514831), 
                          vMr_1 = c(59.9615540623953, 56.4344038234309, 42.3258028675731, 88.1787559741107), 
                          vK3_1 = c(55.0899146106324, 59.0249085113918, 39.3499390075945, 78.6998780151891), 
                          vMr_2 = c(49.1870090082494, 59.7270823671599, 31.6202200767317, 91.3473024438917))%>%
    mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
    arrange(gene)
  
  it("returns a tibble, not a dataframe", {
    expect_true(is_tibble(cpm))
  })
  it("returns exactly 1000 rows and 5 columns", {
    expect_equal(dim(cpm), c(1000, 5))
  })
  it("should contain all the genes found in the original", {
    expect_true(all(td2$gene %in% cpm$gene))
  })
  
  function_cpm <- cpm %>% 
    filter(gene %in% sample_answers$gene) %>%
    arrange(gene)
  
  it("should correctly calculate CPM, testing a subset of known values", {
    expect_true(all(sample_answers == function_cpm))
  })
})

describe("deseq_normalize()", {
  meta_data <- tibble(sample=c("vK3_2", "vMr_1","vK3_1", "vMr_2"), timepoint=c("K3", "Mr", "K3", "Mr"), replicate=c("2","1","1","2"))
  
  normalized <- deseq_normalize(td2, meta_data) %>%
                                   mutate(across(where(is.numeric), \(x) round(x, 3)))
  
  sample_answers <-tibble(gene = c("DYOG91349", "MCAB52159", "QUIS14366", "DRWM24408"), 
                          vK3_2 =  c(23.1379797208499, 13.0779885378717, 25.1499779574455, 456.723599707211), 
                          vMr_1 = c(16.5231110249736, 9.71947707351388, 30.130378927893, 489.861644505099), 
                          vK3_1 = c(17.4371217965184, 9.23141742168622, 25.6428261713506, 413.362357882172), 
                          vMr_2 =  c(15.6092727506343, 15.6092727506343, 40.9743409704149, 458.522387049882)) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
    arrange(gene)
  
  it("returns a tibble, not a dataframe", {
    expect_true(is_tibble(normalized))
  })
  it("returns a tibble with exactly 1000 rows and 5 columns", {
    expect_equal(dim(normalized), c(1000, 5))
  })
  it("should have the same number of genes as the raw data", {
    expect_true(all(td2$gene %in% normalized$gene))
  })
  
  function_deseq <- normalized %>% 
    filter(gene %in% sample_answers$gene) %>%
    arrange(gene)
  
  it("should correctly perform DESeq2 normalization, testing a subset of known values", {
    expect_true(all(function_deseq == sample_answers))
  })
})