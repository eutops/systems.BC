library(rCharts)
# x: buccal
# y: cervical 


# HYPO - HYPO -----------------------------------

hypo_hypo <- dmps %>%
  dplyr::filter(type.y < 0.05 & type.x < 0.05) %>%
  dplyr::filter(dir_type.y < 0 & dir_type.x < 0) %>%
  dplyr::arrange(type.y + type.x) %>%
  dplyr::select(cg) %>%
  dplyr::slice_head(n = 1000)


# Specify the file path where you want to save the .txt file
file_path <- file.path(here("temp"),"hypo_hypo.txt")

# Write the column to a .txt file
write.table(hypo_hypo$cg, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# HYPER - HYPER -----------------------------------

hyper_hyper <- dmps %>%
  dplyr::filter(type.y < 0.05 & type.x < 0.05) %>%
  dplyr::filter(dir_type.y > 0 & dir_type.x > 0) %>%
  dplyr::arrange(type.y + type.x) %>%
  dplyr::select(cg) %>%
  dplyr::slice_head(n = 1000)


# Specify the file path where you want to save the .txt file
file_path <- file.path(here("temp"),"hyper_hyper.txt")

# Write the column to a .txt file
write.table(hyper_hyper$cg, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# HYPER - HYPO -----------------------------------

hyper_hypo <- dmps %>%
  dplyr::filter(type.y < 0.05 & type.x < 0.05) %>%
  dplyr::filter(dir_type.y < 0 & dir_type.x > 0) %>%
  dplyr::arrange(type.y + type.x) %>%
  dplyr::select(cg) %>%
  dplyr::slice_head(n = 1000)


# Specify the file path where you want to save the .txt file
file_path <- file.path(here("temp"),"hyper_hypo.txt")

# Write the column to a .txt file
write.table(hyper_hypo$cg, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# HYPO - HYPER -----------------------------------

hypo_hyper <- dmps %>%
  dplyr::filter(type.y < 0.05 & type.x < 0.05) %>%
  dplyr::filter(dir_type.y > 0 & dir_type.x < 0) %>%
  dplyr::arrange(type.y + type.x) %>%
  dplyr::select(cg) %>%
  dplyr::slice_head(n = 1000)


# Specify the file path where you want to save the .txt file
file_path <- file.path(here("temp"),"hypo_hyper.txt")

# Write the column to a .txt file
write.table(hypo_hyper$cg, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create HTML charts from eforge output files-----------------------------------

# HYPO - HYPO  -----------------------------------


results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hypo.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

color.axis.palette = c();
if (length(which(results$Class2 == 1)) > 0 ) {
  color.axis.palette = c('red');
}
if (length(which(results$Class2 == 2)) > 0 ) {
  color.axis.palette = c(color.axis.palette, '#FF82ab');
}
color.axis.palette = c(color.axis.palette, 'lightblue');
if (length(color.axis.palette) <= 2) {
  color.axis.palette = c(color.axis.palette, 'lightblue'); # Add it twice to force the color if only non-significant values
}

results$log10pvalue <- -log10(results$Pvalue)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Count number of cell types for each tissue (to be able to draw the vertical separation lines afterwards
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

dplot.height=900
dplot.width=2000
bounds.x=60
bounds.y=50
bounds.height=dplot.height - 300
bounds.width=dplot.width - bounds.x - 20

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'log10pvalue',
  x = c('TissueCell'),
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe', 'Class2'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.Unnamed.850k.erc2-DHS'
)

# Force the order on the X-axis
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell',       orderRule = tissue.cell.order[,2])
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'TissueCell', orderRule = as.factor(tissue.cell.order2))

d1$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value

d1$colorAxis(
  type = 'addColorAxis',
  colorSeries = 'Class2',
  palette = color.axis.palette)

# Builds a JS string to add labels for tissues
labels.string = paste(paste0("
    // Adds labels for tissues
    myChart.svg.insert('text', 'g')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-size', 16)
      .attr('font-family', 'Arial')
      .style('fill', '#A9A9A9')
      .attr('transform', 'translate(", (bounds.x + 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), ", 60) rotate(-90)')
      .attr('text-anchor', 'end')
      .text('", names(tissues[2:length(tissues)]), "')
"), collapse='')

# Builds a JS string to add vertical lines to separate tissues
lines.string = paste(paste0("
    // Adds vertical lines between tissues
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y1', 50)
      .attr('x2', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y2', ", (50 + bounds.height), ")
      .style('stroke', '#A9A9A9')
      .style('stroke-dasharray', '10,3,3,3')
"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1$setTemplate(afterScript = paste0("
  <script>
    myChart.draw()
    
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] < 0.01); })
      .style('fill', function(d) {
        var datatype = d['aggField'][4];
        // http://pinetools.com/darken-color to darken hex values for chromatin states
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#ff0000';
        }
      })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#7f0000';
          case 'TssAFlnk':
            return '#7f2200';
          case 'TxFlnk':
            return '#186618';
          case 'Tx':
            return '#004000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#003200';
          case 'Enh':
            return '#727200';
          case 'ZNF-Rpts':
            return '#1e3d32';
          case 'Het':
            return '#31387b';
          case 'TssBiv':
            return '#712222';
          case 'BivFlnk':
            return '#983919';
          case 'EnhBiv':
            return '#66622d';
          case 'ReprPC':
            return '#404040';
          case 'ReprPCWk':
            return '#606060';
          case 'Quies':
            return '#7c7c7c';
          default:
            return '#ff0000';
        }
      });
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return ((d['aggField'][3] >= 0.01) && (d['aggField'][3] < 0.05)); })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#FF82ab';
        }
      })
      .style('stroke-width', '2')
      .style('fill', 'white');
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] >= 0.05); })
      .style('stroke', 'rgb(173, 216, 230)')
      .style('fill', 'white');

    // Substitutes TissueCell labels in X-axis by Cell labels
    myChart.axes[1].shapes
      .selectAll('text')
      .text(function (d) { 
           var i;
           for (i = 0; i < data.length; i += 1) {
               if (data[i].TissueCell === d) {
                   return data[i].Cell;
               }
           }
      })
      .style('text-anchor', 'start');

//       .attr('transform', function() {
//         var t = d3.transform(d3.select(this).attr('transform'));
//         return 'rotate(-90) translate(' + -1*parseFloat(t.translate[0]) + ', ' + -0.7*parseFloat(t.translate[1]) + ')';
//       });

    // Adds title for X-axis
    myChart.axes[1].titleShape
        .style('font-size', 20)

    // Adds title for Y-axis
    myChart.axes[2].titleShape
        .style('font-size', 20)
        .text('-log10 binomial p-value')

    // Adds main title
    myChart.svg.append('text')
      .attr('x', ", (dplot.width / 2), ")
      .attr('y', ", (bounds.y / 2), ")
      .attr('font-size', 24)
      .attr('font-family', 'Arial')
      .attr('font-weight', 'bold')
      .style('fill', 'black')
      .attr('text-anchor', 'middle')
      .text('hypo-hypo')
    ", labels.string, "
    ", lines.string, "
    // Adds vertical line at the far right of the plot
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width), ")
      .attr('y1', ", bounds.y, ")
      .attr('x2', ", (bounds.x + bounds.width), ")
      .attr('y2', ", (bounds.y + bounds.height), ")
      .style('stroke', 'rgb(0,0,0)')
  </script>
"))


d1$save(here("4-markdown-figures","separate-plots","hypo-hypo.html"), cdn = F)

# HYPER - HYPER  -----------------------------------


results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hyper.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

color.axis.palette = c();
if (length(which(results$Class2 == 1)) > 0 ) {
  color.axis.palette = c('red');
}
if (length(which(results$Class2 == 2)) > 0 ) {
  color.axis.palette = c(color.axis.palette, '#FF82ab');
}
color.axis.palette = c(color.axis.palette, 'lightblue');
if (length(color.axis.palette) <= 2) {
  color.axis.palette = c(color.axis.palette, 'lightblue'); # Add it twice to force the color if only non-significant values
}

results$log10pvalue <- -log10(results$Pvalue)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Count number of cell types for each tissue (to be able to draw the vertical separation lines afterwards
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

dplot.height=900
dplot.width=2000
bounds.x=60
bounds.y=50
bounds.height=dplot.height - 300
bounds.width=dplot.width - bounds.x - 20

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'log10pvalue',
  x = c('TissueCell'),
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe', 'Class2'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.Unnamed.850k.erc2-DHS'
)

# Force the order on the X-axis
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell',       orderRule = tissue.cell.order[,2])
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'TissueCell', orderRule = as.factor(tissue.cell.order2))

d1$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value

d1$colorAxis(
  type = 'addColorAxis',
  colorSeries = 'Class2',
  palette = color.axis.palette)

# Builds a JS string to add labels for tissues
labels.string = paste(paste0("
    // Adds labels for tissues
    myChart.svg.insert('text', 'g')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-size', 16)
      .attr('font-family', 'Arial')
      .style('fill', '#A9A9A9')
      .attr('transform', 'translate(", (bounds.x + 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), ", 60) rotate(-90)')
      .attr('text-anchor', 'end')
      .text('", names(tissues[2:length(tissues)]), "')
"), collapse='')

# Builds a JS string to add vertical lines to separate tissues
lines.string = paste(paste0("
    // Adds vertical lines between tissues
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y1', 50)
      .attr('x2', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y2', ", (50 + bounds.height), ")
      .style('stroke', '#A9A9A9')
      .style('stroke-dasharray', '10,3,3,3')
"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1$setTemplate(afterScript = paste0("
  <script>
    myChart.draw()
    
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] < 0.01); })
      .style('fill', function(d) {
        var datatype = d['aggField'][4];
        // http://pinetools.com/darken-color to darken hex values for chromatin states
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#ff0000';
        }
      })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#7f0000';
          case 'TssAFlnk':
            return '#7f2200';
          case 'TxFlnk':
            return '#186618';
          case 'Tx':
            return '#004000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#003200';
          case 'Enh':
            return '#727200';
          case 'ZNF-Rpts':
            return '#1e3d32';
          case 'Het':
            return '#31387b';
          case 'TssBiv':
            return '#712222';
          case 'BivFlnk':
            return '#983919';
          case 'EnhBiv':
            return '#66622d';
          case 'ReprPC':
            return '#404040';
          case 'ReprPCWk':
            return '#606060';
          case 'Quies':
            return '#7c7c7c';
          default:
            return '#ff0000';
        }
      });
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return ((d['aggField'][3] >= 0.01) && (d['aggField'][3] < 0.05)); })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#FF82ab';
        }
      })
      .style('stroke-width', '2')
      .style('fill', 'white');
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] >= 0.05); })
      .style('stroke', 'rgb(173, 216, 230)')
      .style('fill', 'white');

    // Substitutes TissueCell labels in X-axis by Cell labels
    myChart.axes[1].shapes
      .selectAll('text')
      .text(function (d) { 
           var i;
           for (i = 0; i < data.length; i += 1) {
               if (data[i].TissueCell === d) {
                   return data[i].Cell;
               }
           }
      })
      .style('text-anchor', 'start');

//       .attr('transform', function() {
//         var t = d3.transform(d3.select(this).attr('transform'));
//         return 'rotate(-90) translate(' + -1*parseFloat(t.translate[0]) + ', ' + -0.7*parseFloat(t.translate[1]) + ')';
//       });

    // Adds title for X-axis
    myChart.axes[1].titleShape
        .style('font-size', 20)

    // Adds title for Y-axis
    myChart.axes[2].titleShape
        .style('font-size', 20)
        .text('-log10 binomial p-value')

    // Adds main title
    myChart.svg.append('text')
      .attr('x', ", (dplot.width / 2), ")
      .attr('y', ", (bounds.y / 2), ")
      .attr('font-size', 24)
      .attr('font-family', 'Arial')
      .attr('font-weight', 'bold')
      .style('fill', 'black')
      .attr('text-anchor', 'middle')
      .text('hyper-hyper')
    ", labels.string, "
    ", lines.string, "
    // Adds vertical line at the far right of the plot
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width), ")
      .attr('y1', ", bounds.y, ")
      .attr('x2', ", (bounds.x + bounds.width), ")
      .attr('y2', ", (bounds.y + bounds.height), ")
      .style('stroke', 'rgb(0,0,0)')
  </script>
"))


d1$save(here("4-markdown-figures","separate-plots","hyper-hyper.html"), cdn = F)

# HYPO - HYPER  -----------------------------------


results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hyper.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

color.axis.palette = c();
if (length(which(results$Class2 == 1)) > 0 ) {
  color.axis.palette = c('red');
}
if (length(which(results$Class2 == 2)) > 0 ) {
  color.axis.palette = c(color.axis.palette, '#FF82ab');
}
color.axis.palette = c(color.axis.palette, 'lightblue');
if (length(color.axis.palette) <= 2) {
  color.axis.palette = c(color.axis.palette, 'lightblue'); # Add it twice to force the color if only non-significant values
}

results$log10pvalue <- -log10(results$Pvalue)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Count number of cell types for each tissue (to be able to draw the vertical separation lines afterwards
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

dplot.height=900
dplot.width=2000
bounds.x=60
bounds.y=50
bounds.height=dplot.height - 300
bounds.width=dplot.width - bounds.x - 20

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'log10pvalue',
  x = c('TissueCell'),
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe', 'Class2'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.Unnamed.850k.erc2-DHS'
)

# Force the order on the X-axis
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell',       orderRule = tissue.cell.order[,2])
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'TissueCell', orderRule = as.factor(tissue.cell.order2))

d1$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value

d1$colorAxis(
  type = 'addColorAxis',
  colorSeries = 'Class2',
  palette = color.axis.palette)

# Builds a JS string to add labels for tissues
labels.string = paste(paste0("
    // Adds labels for tissues
    myChart.svg.insert('text', 'g')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-size', 16)
      .attr('font-family', 'Arial')
      .style('fill', '#A9A9A9')
      .attr('transform', 'translate(", (bounds.x + 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), ", 60) rotate(-90)')
      .attr('text-anchor', 'end')
      .text('", names(tissues[2:length(tissues)]), "')
"), collapse='')

# Builds a JS string to add vertical lines to separate tissues
lines.string = paste(paste0("
    // Adds vertical lines between tissues
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y1', 50)
      .attr('x2', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y2', ", (50 + bounds.height), ")
      .style('stroke', '#A9A9A9')
      .style('stroke-dasharray', '10,3,3,3')
"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1$setTemplate(afterScript = paste0("
  <script>
    myChart.draw()
    
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] < 0.01); })
      .style('fill', function(d) {
        var datatype = d['aggField'][4];
        // http://pinetools.com/darken-color to darken hex values for chromatin states
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#ff0000';
        }
      })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#7f0000';
          case 'TssAFlnk':
            return '#7f2200';
          case 'TxFlnk':
            return '#186618';
          case 'Tx':
            return '#004000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#003200';
          case 'Enh':
            return '#727200';
          case 'ZNF-Rpts':
            return '#1e3d32';
          case 'Het':
            return '#31387b';
          case 'TssBiv':
            return '#712222';
          case 'BivFlnk':
            return '#983919';
          case 'EnhBiv':
            return '#66622d';
          case 'ReprPC':
            return '#404040';
          case 'ReprPCWk':
            return '#606060';
          case 'Quies':
            return '#7c7c7c';
          default:
            return '#ff0000';
        }
      });
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return ((d['aggField'][3] >= 0.01) && (d['aggField'][3] < 0.05)); })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#FF82ab';
        }
      })
      .style('stroke-width', '2')
      .style('fill', 'white');
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] >= 0.05); })
      .style('stroke', 'rgb(173, 216, 230)')
      .style('fill', 'white');

    // Substitutes TissueCell labels in X-axis by Cell labels
    myChart.axes[1].shapes
      .selectAll('text')
      .text(function (d) { 
           var i;
           for (i = 0; i < data.length; i += 1) {
               if (data[i].TissueCell === d) {
                   return data[i].Cell;
               }
           }
      })
      .style('text-anchor', 'start');

//       .attr('transform', function() {
//         var t = d3.transform(d3.select(this).attr('transform'));
//         return 'rotate(-90) translate(' + -1*parseFloat(t.translate[0]) + ', ' + -0.7*parseFloat(t.translate[1]) + ')';
//       });

    // Adds title for X-axis
    myChart.axes[1].titleShape
        .style('font-size', 20)

    // Adds title for Y-axis
    myChart.axes[2].titleShape
        .style('font-size', 20)
        .text('-log10 binomial p-value')

    // Adds main title
    myChart.svg.append('text')
      .attr('x', ", (dplot.width / 2), ")
      .attr('y', ", (bounds.y / 2), ")
      .attr('font-size', 24)
      .attr('font-family', 'Arial')
      .attr('font-weight', 'bold')
      .style('fill', 'black')
      .attr('text-anchor', 'middle')
      .text('hypo-hyper')
    ", labels.string, "
    ", lines.string, "
    // Adds vertical line at the far right of the plot
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width), ")
      .attr('y1', ", bounds.y, ")
      .attr('x2', ", (bounds.x + bounds.width), ")
      .attr('y2', ", (bounds.y + bounds.height), ")
      .style('stroke', 'rgb(0,0,0)')
  </script>
"))


d1$save(here("4-markdown-figures","separate-plots","hypo-hyper.html"), cdn = F)

# HYPER - HYPO  -----------------------------------


results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hypo.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

color.axis.palette = c();
if (length(which(results$Class2 == 1)) > 0 ) {
  color.axis.palette = c('red');
}
if (length(which(results$Class2 == 2)) > 0 ) {
  color.axis.palette = c(color.axis.palette, '#FF82ab');
}
color.axis.palette = c(color.axis.palette, 'lightblue');
if (length(color.axis.palette) <= 2) {
  color.axis.palette = c(color.axis.palette, 'lightblue'); # Add it twice to force the color if only non-significant values
}

results$log10pvalue <- -log10(results$Pvalue)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Count number of cell types for each tissue (to be able to draw the vertical separation lines afterwards
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

dplot.height=900
dplot.width=2000
bounds.x=60
bounds.y=50
bounds.height=dplot.height - 300
bounds.width=dplot.width - bounds.x - 20

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'log10pvalue',
  x = c('TissueCell'),
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe', 'Class2'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.Unnamed.850k.erc2-DHS'
)

# Force the order on the X-axis
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell',       orderRule = tissue.cell.order[,2])
d1$xAxis( type = 'addCategoryAxis', grouporderRule = 'TissueCell', orderRule = as.factor(tissue.cell.order2))

d1$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value

d1$colorAxis(
  type = 'addColorAxis',
  colorSeries = 'Class2',
  palette = color.axis.palette)

# Builds a JS string to add labels for tissues
labels.string = paste(paste0("
    // Adds labels for tissues
    myChart.svg.insert('text', 'g')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-size', 16)
      .attr('font-family', 'Arial')
      .style('fill', '#A9A9A9')
      .attr('transform', 'translate(", (bounds.x + 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), ", 60) rotate(-90)')
      .attr('text-anchor', 'end')
      .text('", names(tissues[2:length(tissues)]), "')
"), collapse='')

# Builds a JS string to add vertical lines to separate tissues
lines.string = paste(paste0("
    // Adds vertical lines between tissues
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y1', 50)
      .attr('x2', ", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), ")
      .attr('y2', ", (50 + bounds.height), ")
      .style('stroke', '#A9A9A9')
      .style('stroke-dasharray', '10,3,3,3')
"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1$setTemplate(afterScript = paste0("
  <script>
    myChart.draw()
    
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] < 0.01); })
      .style('fill', function(d) {
        var datatype = d['aggField'][4];
        // http://pinetools.com/darken-color to darken hex values for chromatin states
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#ff0000';
        }
      })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#7f0000';
          case 'TssAFlnk':
            return '#7f2200';
          case 'TxFlnk':
            return '#186618';
          case 'Tx':
            return '#004000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#003200';
          case 'Enh':
            return '#727200';
          case 'ZNF-Rpts':
            return '#1e3d32';
          case 'Het':
            return '#31387b';
          case 'TssBiv':
            return '#712222';
          case 'BivFlnk':
            return '#983919';
          case 'EnhBiv':
            return '#66622d';
          case 'ReprPC':
            return '#404040';
          case 'ReprPCWk':
            return '#606060';
          case 'Quies':
            return '#7c7c7c';
          default:
            return '#ff0000';
        }
      });
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return ((d['aggField'][3] >= 0.01) && (d['aggField'][3] < 0.05)); })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#FF82ab';
        }
      })
      .style('stroke-width', '2')
      .style('fill', 'white');
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] >= 0.05); })
      .style('stroke', 'rgb(173, 216, 230)')
      .style('fill', 'white');

    // Substitutes TissueCell labels in X-axis by Cell labels
    myChart.axes[1].shapes
      .selectAll('text')
      .text(function (d) { 
           var i;
           for (i = 0; i < data.length; i += 1) {
               if (data[i].TissueCell === d) {
                   return data[i].Cell;
               }
           }
      })
      .style('text-anchor', 'start');

//       .attr('transform', function() {
//         var t = d3.transform(d3.select(this).attr('transform'));
//         return 'rotate(-90) translate(' + -1*parseFloat(t.translate[0]) + ', ' + -0.7*parseFloat(t.translate[1]) + ')';
//       });

    // Adds title for X-axis
    myChart.axes[1].titleShape
        .style('font-size', 20)

    // Adds title for Y-axis
    myChart.axes[2].titleShape
        .style('font-size', 20)
        .text('-log10 binomial p-value')

    // Adds main title
    myChart.svg.append('text')
      .attr('x', ", (dplot.width / 2), ")
      .attr('y', ", (bounds.y / 2), ")
      .attr('font-size', 24)
      .attr('font-family', 'Arial')
      .attr('font-weight', 'bold')
      .style('fill', 'black')
      .attr('text-anchor', 'middle')
      .text('hyper-hypo')
    ", labels.string, "
    ", lines.string, "
    // Adds vertical line at the far right of the plot
    myChart.svg.append('line')
      .attr('x1', ", (bounds.x + bounds.width), ")
      .attr('y1', ", bounds.y, ")
      .attr('x2', ", (bounds.x + bounds.width), ")
      .attr('y2', ", (bounds.y + bounds.height), ")
      .style('stroke', 'rgb(0,0,0)')
  </script>
"))


d1$save(here("4-markdown-figures","separate-plots","hyper-hypo.html"), cdn = F)



# Create PDF charts from eforge output files-----------------------------------

# HYPER - HYPO PDF -----------------------------------

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hypo.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Plot an empty chart first
pdf('hyper_hypo.pdf', width=22.4, height=8)
ymax = max(-log10(results$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5,4,3,1)+0.1)
#plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for erc2-DHS Unnamed',
#    ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)
plot(NA,ylab='', xlab='', main='hyper_hypo',
     ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='grey')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='grey', adj=1, srt=90, cex=1.2) 

# Add points (internal color first)
palette(c('red', 'palevioletred1', 'white'))
points(results$TissueCell, -log10(results$Pvalue), pch=19, col = results$Class2, xaxt='n')

# Add contour to the points
palette(c('black', 'palevioletred1', 'steelblue3'))
points(results$TissueCell, -log10(results$Pvalue), pch=1, col = results$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

# Add legend (internal color first)
palette(c('white', 'palevioletred1', 'red'))
legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('steelblue3', 'palevioletred1', 'black'))
legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

palette('default')
dev.off()





# HYPER - HYPER PDF -----------------------------------

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hyper.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Plot an empty chart first
pdf('hyper_hyper.pdf', width=22.4, height=8)
ymax = max(-log10(results$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5,4,3,1)+0.1)
#plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for erc2-DHS Unnamed',
#    ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)
plot(NA,ylab='', xlab='', main='hyper_hyper',
     ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='grey')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='grey', adj=1, srt=90, cex=1.2) 

# Add points (internal color first)
palette(c('red', 'palevioletred1', 'white'))
points(results$TissueCell, -log10(results$Pvalue), pch=19, col = results$Class2, xaxt='n')

# Add contour to the points
palette(c('black', 'palevioletred1', 'steelblue3'))
points(results$TissueCell, -log10(results$Pvalue), pch=1, col = results$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

# Add legend (internal color first)
palette(c('white', 'palevioletred1', 'red'))
legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('steelblue3', 'palevioletred1', 'black'))
legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

palette('default')
dev.off()




# HYPO - HYPO PDF -----------------------------------

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hypo.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Plot an empty chart first
pdf('hypo_hypo.pdf', width=22.4, height=8)
ymax = max(-log10(results$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5,4,3,1)+0.1)
#plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for erc2-DHS Unnamed',
#    ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)
plot(NA,ylab='', xlab='', main='hypo_hypo',
     ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='grey')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='grey', adj=1, srt=90, cex=1.2) 

# Add points (internal color first)
palette(c('red', 'palevioletred1', 'white'))
points(results$TissueCell, -log10(results$Pvalue), pch=19, col = results$Class2, xaxt='n')

# Add contour to the points
palette(c('black', 'palevioletred1', 'steelblue3'))
points(results$TissueCell, -log10(results$Pvalue), pch=1, col = results$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

# Add legend (internal color first)
palette(c('white', 'palevioletred1', 'red'))
legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('steelblue3', 'palevioletred1', 'black'))
legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

palette('default')
dev.off()



# HYPO - HYPER PDF -----------------------------------

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hyper.tsv.gz")), header = TRUE, sep="\t")

# Class splits the data into non-significant, marginally significant and significant according to 0.05 and 0.01 (in -log10 scale)
results$Class <- cut(results$Pvalue, breaks =c(0, 0.01, 0.05, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results$Class2 <- cut(results$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results$TissueCell <- factor(results$TissueCell, levels=tissue.cell.order2)

# Plot an empty chart first
pdf('hypo_hyper.pdf', width=22.4, height=8)
ymax = max(-log10(results$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5,4,3,1)+0.1)
#plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for erc2-DHS Unnamed',
#    ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)
plot(NA,ylab='', xlab='', main='hypo_hyper',
     ylim=c(ymin,ymax), las=2, pch=19, col = results$Class2, xaxt='n', xlim=c(0,length(levels(results$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='grey')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='grey', adj=1, srt=90, cex=1.2) 

# Add points (internal color first)
palette(c('red', 'palevioletred1', 'white'))
points(results$TissueCell, -log10(results$Pvalue), pch=19, col = results$Class2, xaxt='n')

# Add contour to the points
palette(c('black', 'palevioletred1', 'steelblue3'))
points(results$TissueCell, -log10(results$Pvalue), pch=1, col = results$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

# Add legend (internal color first)
palette(c('white', 'palevioletred1', 'red'))
legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('steelblue3', 'palevioletred1', 'black'))
legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

palette('default')
dev.off()

# TABLES - ALL -----------------------------------


library("writexl")

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hypo.tsv.gz")), header = TRUE, sep="\t")
save(results,file=here("5-supplements","eFORGE","hyper_hypo.Rdata"))
write_xlsx(results, here("5-supplements","eFORGE","hyper_hypo.xlsx"))

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hyper_hyper.tsv.gz")), header = TRUE, sep="\t")
save(results,file=here("5-supplements","eFORGE","hyper_hyper.Rdata"))
write_xlsx(results, here("5-supplements","eFORGE","hyper_hyper.xlsx"))

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hypo.tsv.gz")), header = TRUE, sep="\t")
save(results,file=here("5-supplements","eFORGE","hypo_hypo.Rdata"))
write_xlsx(results, here("5-supplements","eFORGE","hypo_hypo.xlsx"))

results<-read.table(file.path(here("0-data","dataframes-plotting","eforge_hypo_hyper.tsv.gz")), header = TRUE, sep="\t")
save(results,file=here("5-supplements","eFORGE","hypo_hyper.Rdata"))
write_xlsx(results, here("5-supplements","eFORGE","hypo_hyper.xlsx"))



