var dataset;
var toMake = "DN1";


var totalGenes;
const Turbocolors = ["#23171b","#271a28","#2b1c33","#2f1e3f","#32204a","#362354","#39255f","#3b2768","#3e2a72","#402c7b","#422f83","#44318b","#453493","#46369b","#4839a2","#493ca8","#493eaf","#4a41b5","#4a44bb","#4b46c0","#4b49c5","#4b4cca","#4b4ecf","#4b51d3","#4a54d7","#4a56db","#4959de","#495ce2","#485fe5","#4761e7","#4664ea","#4567ec","#446aee","#446df0","#426ff2","#4172f3","#4075f5","#3f78f6","#3e7af7","#3d7df7","#3c80f8","#3a83f9","#3985f9","#3888f9","#378bf9","#368df9","#3590f8","#3393f8","#3295f7","#3198f7","#309bf6","#2f9df5","#2ea0f4","#2da2f3","#2ca5f1","#2ba7f0","#2aaaef","#2aaced","#29afec","#28b1ea","#28b4e8","#27b6e6","#27b8e5","#26bbe3","#26bde1","#26bfdf","#25c1dc","#25c3da","#25c6d8","#25c8d6","#25cad3","#25ccd1","#25cecf","#26d0cc","#26d2ca","#26d4c8","#27d6c5","#27d8c3","#28d9c0","#29dbbe","#29ddbb","#2adfb8","#2be0b6","#2ce2b3","#2de3b1","#2ee5ae","#30e6ac","#31e8a9","#32e9a6","#34eba4","#35eca1","#37ed9f","#39ef9c","#3af09a","#3cf197","#3ef295","#40f392","#42f490","#44f58d","#46f68b","#48f788","#4af786","#4df884","#4ff981","#51fa7f","#54fa7d","#56fb7a","#59fb78","#5cfc76","#5efc74","#61fd71","#64fd6f","#66fd6d","#69fd6b","#6cfd69","#6ffe67","#72fe65","#75fe63","#78fe61","#7bfe5f","#7efd5d","#81fd5c","#84fd5a","#87fd58","#8afc56","#8dfc55","#90fb53","#93fb51","#96fa50","#99fa4e","#9cf94d","#9ff84b","#a2f84a","#a6f748","#a9f647","#acf546","#aff444","#b2f343","#b5f242","#b8f141","#bbf03f","#beef3e","#c1ed3d","#c3ec3c","#c6eb3b","#c9e93a","#cce839","#cfe738","#d1e537","#d4e336","#d7e235","#d9e034","#dcdf33","#dedd32","#e0db32","#e3d931","#e5d730","#e7d52f","#e9d42f","#ecd22e","#eed02d","#f0ce2c","#f1cb2c","#f3c92b","#f5c72b","#f7c52a","#f8c329","#fac029","#fbbe28","#fdbc28","#feb927","#ffb727","#ffb526","#ffb226","#ffb025","#ffad25","#ffab24","#ffa824","#ffa623","#ffa323","#ffa022","#ff9e22","#ff9b21","#ff9921","#ff9621","#ff9320","#ff9020","#ff8e1f","#ff8b1f","#ff881e","#ff851e","#ff831d","#ff801d","#ff7d1d","#ff7a1c","#ff781c","#ff751b","#ff721b","#ff6f1a","#fd6c1a","#fc6a19","#fa6719","#f96418","#f76118","#f65f18","#f45c17","#f25916","#f05716","#ee5415","#ec5115","#ea4f14","#e84c14","#e64913","#e44713","#e24412","#df4212","#dd3f11","#da3d10","#d83a10","#d5380f","#d3360f","#d0330e","#ce310d","#cb2f0d","#c92d0c","#c62a0b","#c3280b","#c1260a","#be2409","#bb2309","#b92108","#b61f07","#b41d07","#b11b06","#af1a05","#ac1805","#aa1704","#a81604","#a51403","#a31302","#a11202","#9f1101","#9d1000","#9b0f00","#9a0e00","#980e00","#960d00","#950c00","#940c00","#930c00","#920c00","#910b00","#910c00","#900c00","#900c00","#900c00"]
var colorMap;


var slideWidth = 200; // default window width for chromosome view

var orthologs = window.orthologs;   //trueMatch made by the collinearity script

window.onresize = function () { location.reload() }



d3.text('canola.gff').
then(function(d){
    var result = "chromosome\tgene\tstart\tend\n" + d;
    dataset = d3.tsvParse(result)
    totalGenes  = dataset.length;
    generateVisualization();

}) 


var windWidth =  window.innerWidth 



function generateVisualization() {
   

    d3.select('.genomeView').remove();
    d3.select('.genomeView2').remove();
    
    var slider = document.getElementById("myRange");
    var output = document.getElementById("range")
    output.innerHTML = slider.value;

    slider.oninput = function() {
        output.innerHTML = this.value;
        
        
      }
    slider.onchange = function(){
        d3.select('.genomeView').remove();
        generateVisualization();
    }
    

    
    var currentChromosome = dataset[0].chromosome


    var chromosomeSorted ={};  //{chromosome: {geneName: [geneStart: geneEnd]....}....}


    for (let entry of dataset){
        var gene = entry.gene;

        if (chromosomeSorted[entry.chromosome]){
        chromosomeSorted[entry.chromosome][gene] = [entry.start, entry.end]}
        else{
            chromosomeSorted[entry.chromosome] =  {}
            chromosomeSorted[entry.chromosome][gene] = [entry.start, entry.end]
        }

    }


    var chromDrawarray = makePopulationData(chromosomeSorted, slider.value); //Array of linear data to draw


    //Total lengths for scaling
    var combinedChromosomeLengthDN = 0;
    var combinedChromosomeLengthN = 0;


    chromosomeLengths = {};
    colorMap = {};
    var counter = 1;



    var cNum =  0
    for ( let chromosome of Object.keys(chromDrawarray)){
        counter++
        colorMap[chromosome] = Turbocolors[counter*Math.round(Turbocolors.length/40)];
        var currentChromLength = getChromosomeLength(chromosome,chromosomeSorted);
        chromosomeLengths[chromosome] = currentChromLength;

        if(chromosome.slice(0,1)=='D'){
            combinedChromosomeLengthDN  +=currentChromLength;
            }else{
            combinedChromosomeLengthN += currentChromLength;
            }

    }

    //GenomeView row 1
    d3.select('#wholeScreen').append('div').attr("class","genomeView text-center").attr("style","position: relative;");
    for ( let chromosome of Object.keys(chromDrawarray).slice(0,19).sort(sortAlphaNum)){
        
        cNum++;

        numGenesArray = chromDrawarray[chromosome]
        

        
        var widthScale = d3.scaleLinear()
        .domain([0,combinedChromosomeLengthDN])
    
        .range([0, windWidth-100]);


        var w = widthScale(chromosomeLengths[chromosome]),
            h = '100';
            

        var max = Math.max.apply(null, numGenesArray);
        var min = Math.min.apply(null, numGenesArray);
        
        var colorScale = d3.scaleLinear()
        .domain([min,max])
    
        .range(["white", colorMap[chromosome]]);



        var svg = d3.select(".genomeView").append("svg").attr("width", w).attr("height", h).attr("id",chromosome).on("click",function(){

            makeBigChrom(this.id, chromDrawarray, d3.scaleLinear()
            .domain([min,max])
        
            .range(["white", colorMap[this.id]]),chromosomeSorted, slider);
            toMake = this.id
        })
        
        


        svg.chrom = chromosome;
        d3.select(".genomeView").append("text")
        .text("      ");
        var dataChrom = chromosomeSorted[chromosome]
        var m = [].concat.apply([], dataChrom);


        var upper = (Math.max.apply(null, m));
        var lower = (Math.min.apply(null, m));
        makeMap(numGenesArray, w, chromosome,svg,colorScale, 50, false, chromosomeSorted, lower, upper);
            
    }

    //Genome view row 2
    d3.select('#wholeScreen').append('div').attr("class","genomeView2 text-center").attr("style","position: relative;");
    for ( let chromosome of Object.keys(chromDrawarray).slice(19).sort(sortAlphaNum)){
        
        cNum++;

        numGenesArray = chromDrawarray[chromosome]
        

        
        var widthScale = d3.scaleLinear()
        .domain([0,combinedChromosomeLengthN])
    
        .range([0, windWidth-100]);


        var w = widthScale(chromosomeLengths[chromosome]),
            h = '100';
            

        var max = Math.max.apply(null, numGenesArray);
        var min = Math.min.apply(null, numGenesArray);
        
        var colorScale = d3.scaleLinear()
        .domain([min,max])
    
        .range(["white", colorMap[chromosome]]);



        var svg = d3.select(".genomeView2").append("svg").attr("width", w).attr("height", h).attr("id",chromosome).on("click",function(){

            makeBigChrom(this.id, chromDrawarray, d3.scaleLinear()
            .domain([min,max])
        
            .range(["white", colorMap[this.id]]),chromosomeSorted, slider);
            toMake = this.id
        })
        svg.chrom = chromosome;
        d3.select(".genomeView2").append("text")
        .text("      ");
        var dataChrom = chromosomeSorted[chromosome]
        var m = [].concat.apply([], dataChrom);


        var upper = (Math.max.apply(null, m));
        var lower = (Math.min.apply(null, m));
        makeMap(numGenesArray, w, chromosome,svg,colorScale, 50, false, chromosomeSorted, lower, upper);
            
    }

    for ( let chromosome of Object.keys(chromDrawarray)){

    var currentSVG = document.getElementById(chromosome)
    var chromeStart = currentSVG.getBoundingClientRect().left;
    var chromeWidth  = currentSVG.getBoundingClientRect().width;

    if(chromosome.slice(0,1)=='D'){
    d3.select(`.genomeView`).append('div').attr("style",`position: absolute; top: 60px;left: ${chromeStart}; width: ${chromeWidth};`).attr("class", "text-justify").text(chromosome)}
    else{
        d3.select(`.genomeView2`).append('div').attr("style",`position: absolute; top: 60px;left: ${chromeStart}; width: ${chromeWidth};`).attr("class", "text-justify").text(chromosome)
    }

    }






    makeBigChrom(toMake, chromDrawarray, d3.scaleLinear()
    .domain([min,max])

    .range(["white", colorMap[toMake]]), chromosomeSorted, slider)
    var geneSearched = document.getElementById("geneToSearch");
    var button = document.getElementById("geneSearchButton")
    button.onclick = function(){
        

        searchGene(geneSearched.value, chromDrawarray, colorScale, chromosomeSorted, slider)
    }

    
    

}
    

    


function makeBigChrom(chromosome, chromDrawarray, colorScale, chromosomeSorted, slider){
    d3.select('.chromosomePortion').remove();
    d3.select('.chromosomeShower').remove();

    d3.select('#wholeScreen').append('div').attr("class","chromosomeShower text-center").attr("style","position: relative").text(chromosome);
    var svg2 = d3.select(".chromosomeShower").append("svg").attr("width", windWidth).attr("height", 150);
    d3.select(".chromosomeShower").append('div').attr("class", "slider").attr("style",style=`position: absolute; width: ${slideWidth} ; height: 100px; top: 13px; opacity: 0.75; border: 15px solid grey;`).attr("id","slidingWindow")
    

    
    var dataChrom = Object.values(chromosomeSorted[chromosome])
    var m = [].concat.apply([], dataChrom);

        
    var upper = (Math.max.apply(null, m));
    var lower = (Math.min.apply(null, m));

    var max = Math.max.apply(null, chromDrawarray[chromosome]);
    var min = Math.min.apply(null, chromDrawarray[chromosome]);



        
    colorScale = d3.scaleLinear()
    .domain([min,max])

    .range(["white", colorMap[chromosome]]);

    makeMap(chromDrawarray[chromosome], windWidth, chromosome,svg2,colorScale,75, true, chromosomeSorted, lower, upper);

    makeChromPortion(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(0, d3.scaleLinear()
            .domain([0, (chromDrawarray[toMake]).length - 1])
            .range([0, windWidth])), slider)


    const position = { x: 0, y: 0 }

    interact('.slider').draggable({
    listeners: {

        move (event) {
        position.x += event.dx

        event.target.style.transform =
            `translate(${position.x}px, ${position.y}px)`
        },
        end (event){
            makeChromPortion(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(position.x, d3.scaleLinear()
            .domain([0, (chromDrawarray[chromosome]).length - 1])
            .range([0, windWidth])), slider)

        }
    }
    })

    interact('.slider')
  .resizable({
    edges: { top: false, left: false, bottom: false, right: true },
    listeners: {
      move: function (event) {
        let { x, y } = event.target.dataset

        x = (parseFloat(x) || 0) + event.deltaRect.left


        Object.assign(event.target.style, {
          width: `${event.rect.width}px`,
          transform: `translate(${x}px, ${y}px)`
        })

        Object.assign(event.target.dataset, { x, y })
      },
      end (event){
          makeChromPortion(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(position.x, d3.scaleLinear()
          .domain([0, (chromDrawarray[chromosome]).length - 1])
          .range([0, windWidth])), slider)

      }
    }
  })

    


}

function makeMap(n, wid, chromo,svg,colorScale, h, axisRequried, chromosomeSorted, lower, upper){

    var data_length = n.length;
    var xScale = d3.scaleLinear()
            .domain([0, data_length - 1])
            .range([0, wid]);

    var Scale = d3.scaleLinear()
            .domain([lower, upper])
            .range([0, wid]);

    svg.selectAll("rect").data(n).enter().append("rect").on("click",function(d,i){

    })
            .attr("x", function(d, i) {
                return xScale(i);
            }).attr("y", function(d) {
                return 0;
            }) //move away from top 	
            .attr("width", wid/ n.length)
            .attr("height", function(d) {
                return h;
            })
            .attr("fill", function(d) {

                return colorScale(d);
                
            })
            .append("svg:title")
            .text(chromo)
            
    
    if (axisRequried){

    
    

        var x_axis = d3.axisBottom()
    .scale(Scale);

    svg.append("g").attr("transform", "translate(10, 90)")
    .call(x_axis);
    }
    

}


function makeChromPortion(chromosome, chromDrawarray, colorScale, chromosomeSorted, index, slider){


    d3.select('.chromosomePortion').remove();

    d3.select('#wholeScreen').append('div').attr("class","chromosomePortion text-center").attr("style","position: relative");
    d3.select(".chromosomePortion").append('div').attr("class", "slider2").attr("id","geneSelector").attr("style","position: absolute;width: 0px;height: 0px;top: 75px;cursor: pointer;opacity: 0.75;left: 0x; border-bottom: 16px solid black; border-left: 10px solid transparent; border-right: 10px solid transparent; ")
    
    var svg2 = d3.select(".chromosomePortion").append("svg").attr("width", windWidth).attr("height", 150);


    var element = document.getElementById("slidingWindow")
    slideWidth = parseInt(element.style.width)
        
    var lower = (index*slider.value);
    var right = sliderToIndex(slideWidth, d3.scaleLinear()
            .domain([0, (chromDrawarray[chromosome]).length - 1])
            .range([0, windWidth]))
    var upper = ((index+right)*slider.value);
    var requiredGenes = {};
    requiredGenes[chromosome] = [];
    for ( let gene of Object.values(chromosomeSorted[chromosome])){
        if (gene[0] >= lower && gene[1] <= upper){

            requiredGenes[chromosome].push(gene)
        }
    }
    var m = [].concat.apply([], requiredGenes[chromosome]);
    var max = Math.max.apply(null, m);
    var min = Math.min.apply(null, m);

    var upperLimit = String(Math.max.apply(null, m));
    var currentBase = String(Math.min.apply(null, m));



    var num  = Math.floor((upperLimit - currentBase)/5*windWidth)

    var PopulationData = makePopulationData(requiredGenes, 5*windWidth)

    var reqData = PopulationData[chromosome]

    



    var merged = [].concat.apply([], reqData);
    var max = Math.max.apply(null, merged);
    var min = Math.min.apply(null, merged);
    
    
    var colorScale = d3.scaleLinear()
    .domain([min,max])

    .range(["white", "orange"]);
    
    makeMap(reqData, windWidth, chromosome,svg2,colorScale,75, true, chromosomeSorted, lower, upper);
    makeGenemap(chromosome, chromDrawarray, colorScale, chromosomeSorted, 0,lower, upper, slider)

    const position = { x: 0, y: 0 }

    interact('.slider2').draggable({
    listeners: {

        move (event) {
        position.x += event.dx

        event.target.style.transform =
            `translate(${position.x}px, ${position.y}px)`


            var element = document.getElementById("slidingWindow")
            var slideWidth = parseInt(element.style.width)
            makeGenemap(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(position.x , d3.scaleLinear()
            .domain([0, slideWidth])
            .range([0, windWidth])),lower,upper, slider)
        },
        end (event){
            

        }
    }
    })
    

}

function makeGenemap(chromosome, chromDrawarray, colorScale, chromosomeSorted, index,start, end, slider){
    d3.select('.geneMap').remove();

    d3.select('#wholeScreen').append('div').attr("class","geneMap text-center").attr("style","position: relative;");

    
    var svg2 = d3.select(".geneMap").append("svg").attr("width", windWidth).attr("height", 150);

    var element = document.getElementById("slidingWindow")
    var slideWidth = parseInt(element.style.width)

    var point = sliderToIndex(index, d3.scaleLinear()
    .domain([start, end])
    .range([0, slideWidth]));

    
    
    var startPoint = Math.max(0, point - 10000);
    var endPoint =  point + 10000;



    var toDraw =  [];
    var corresPondingGene = [];
    for ( let i of Object.keys(chromosomeSorted[chromosome])){
        if (chromosomeSorted[chromosome][i][0] > startPoint && chromosomeSorted[chromosome][i][1]<endPoint){
            toDraw.push(chromosomeSorted[chromosome][i])
            corresPondingGene.push(i)
        }
    }


    var xScale = d3.scaleLinear()
            .domain([startPoint, endPoint])
            .range([0, windWidth]);

    var Scale = d3.scaleLinear()
            .domain([0, endPoint -startPoint])
            .range([0, windWidth]);

    svg2.selectAll("rect").data(toDraw).enter().append("rect")
            .attr("x", function(d, i) {
                return xScale(d[0]);
            }).attr("y", function(d) {
                return 0;
            }) //move away from top 	
            .attr("width", function(d,i){

                return Scale(d[1]-d[0])
            })
            .attr("height", function(d) {
                return 75;
            })
            .attr("fill", function(d) {

                return "#006666";
                
            })
            .attr("id",function(d,i) {

                return corresPondingGene[i];;
                
            })
            
            .on("click",function(d,i){


                searchGene(this.id, chromDrawarray, colorScale, chromosomeSorted, slider);

            })

            .append("svg:title")
            .text(function(d, i) { 

                return `${corresPondingGene[i]}`; })
            

    var x_axis = d3.axisBottom()
    .scale(xScale);

    svg2.append("g").attr("transform", "translate(0, 100)")
    .call(x_axis);



}



function sliderToIndex(x, chartScale) {
    
    const location = Math.abs(x)

    return Math.round(chartScale.invert(location))

}

function searchGene(geneSearched, chromDrawarray, colorScale, chromosomeSorted, slider){



    d3.select('#orthologToggler').remove();

    var allOrthologs = getOrthologs(geneSearched)

    makeMarkers(allOrthologs, chromosomeSorted, geneSearched)
    var counter  = 0;
    
    if (allOrthologs.length>1){

        var element = document.getElementById("homolog-panel");
        element.innerHTML += '<div id="orthologToggler"></div>'
        var el = document.getElementById("orthologToggler")
        // (document.getElementById("previousOrtholog")).remove();
        // element.removeChild(document.getElementById("nextOrtholog"));
        el.innerHTML += '<button type="button" id="previousOrtholog" class="btn btn-outline-primary">Previous Ortholog</button>'
        el.innerHTML += '<button type="button" id="nextOrtholog" class="btn btn-outline-primary">Next Ortholog</button>'
    }

    focusedGene(geneSearched, chromDrawarray, colorScale, chromosomeSorted, slider);

    var moveTonextOrtholog = document.getElementById("nextOrtholog")
    if (moveTonextOrtholog){
    moveTonextOrtholog.onclick = function(){
        counter++;
        if (counter < allOrthologs.length){
            focusedGene(allOrthologs[counter], chromDrawarray, colorScale, chromosomeSorted, slider);
        }else{
            counter = allOrthologs.length-1
        }
    }
    }

    
    var moveTopreviousOrtholog = document.getElementById("previousOrtholog")
    if (moveTopreviousOrtholog){
    moveTopreviousOrtholog.onclick = function(){
        counter--;
        if (counter >= 0){
            focusedGene(allOrthologs[counter], chromDrawarray, colorScale, chromosomeSorted, slider);
            
        }else{
            counter = 0
        }
    }
    }
}

function  findGene(geneSearched) {


    for ( let entry of dataset){
        if (entry.gene.toLowerCase() == geneSearched.toLowerCase()){

            return entry;
        }
    }
}

function focusedGene(geneSearched, chromDrawarray, colorScale, chromosomeSorted, slider){
    var gene = findGene(geneSearched);
    var chromosome = gene.chromosome;
    var start = gene.start;
    var end = gene.end;

    makeBigChrom(chromosome, chromDrawarray, colorScale, chromosomeSorted, slider)
    var element = document.getElementById("slidingWindow")
    var slideWidth = parseInt(element.style.width)


    var data = chromDrawarray[chromosome];
    var data_length = data.length;

    var xScale = d3.scaleLinear()
            .domain([0, data_length - 1])
            .range([0, windWidth]);


    var slidingWindowPosition = xScale(start/slider.value) - (slideWidth/2);
    if (slidingWindowPosition<0){
        slidingWindowPosition = 0
    }

    element.style.transform =
            `translate(${slidingWindowPosition}px, 0px)`
    makeChromPortion(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(xScale(start/slider.value) - (slideWidth/2), d3.scaleLinear()
            .domain([0, (chromDrawarray[chromosome]).length - 1])
            .range([0, windWidth])), slider)

    
    

    var index = sliderToIndex(slidingWindowPosition, d3.scaleLinear()
    .domain([0, (chromDrawarray[chromosome]).length - 1])
    .range([0, windWidth]))
    var lower = (index*slider.value);


    var element = document.getElementById("slidingWindow")
    var slideWidth = parseInt(element.style.width)

    var right = sliderToIndex(slideWidth, xScale)
    var upper = ((index+right)*slider.value);

    var Scale2 = d3.scaleLinear()
            .domain([lower, upper])
            .range([0, windWidth]);

    var element2 = document.getElementById("geneSelector")


    toGo = Scale2(start)

    element2.style.transform =
            `translate(${toGo}px, 0px)`

    makeGenemapConditional(chromosome, chromDrawarray, colorScale, chromosomeSorted, sliderToIndex(toGo, d3.scaleLinear()
    .domain([0, slideWidth])
    .range([0, windWidth])),lower, upper, slider, start, end, geneSearched )
    

}

function makeGenemapConditional(chromosome, chromDrawarray, colorScale, chromosomeSorted, index,start, end, slider, geneStart, geneEnd, geneSearched){
    d3.select('.geneMap').remove();

    d3.select('#wholeScreen').append('div').attr("class","geneMap text-center").attr("style","position: relative;")

    
    var svg2 = d3.select(".geneMap").append("svg").attr("width", windWidth).attr("height", 150);

    var element = document.getElementById("slidingWindow")
    var slideWidth = parseInt(element.style.width)

    var point = sliderToIndex(index, d3.scaleLinear()
    .domain([start, end])
    .range([0, slideWidth]));

    
    
    var startPoint = Math.max(0, point - 10000);
    var endPoint =  point + 10000;



    var toDraw =  [];
    var corresPondingGene = [];

    for (let i of Object.keys(chromosomeSorted[chromosome])){
        if (chromosomeSorted[chromosome][i][0] > startPoint && chromosomeSorted[chromosome][i][1]<endPoint){
            toDraw.push(chromosomeSorted[chromosome][i])
            corresPondingGene.push(i)
        }
    }


    var xScale = d3.scaleLinear()
            .domain([startPoint, endPoint])
            .range([0, windWidth]);

    var Scale = d3.scaleLinear()
            .domain([0, endPoint -startPoint])
            .range([0, windWidth]);

    svg2.selectAll("rect").data(toDraw).enter().append("rect")
            .attr("x", function(d, i) {
                return xScale(d[0]);
            }).attr("y", function(d) {
                return 0;
            }) //move away from top 	
            .attr("width", function(d,i){

                return Scale(d[1]-d[0])
            })
            .attr("height", function(d) {
                return 75;
            })
            .attr("fill", function(d) {

                if (d[1] == geneEnd && d[0] == geneStart){
                    return "#d42447";
                }
                else{
                     return "#006666";
                }
            }).attr("id",function(d,i) {

                return corresPondingGene[i];;
                
            })
            
            .on("click",function(d,i){


                searchGene(this.id, chromDrawarray, colorScale, chromosomeSorted, slider);

            })

            .append("svg:title")
            .text(function(d, i) { return `${corresPondingGene[i]} \n Orthologs:  ${getOrthologs(corresPondingGene[i]).length}`; })

            

    var x_axis = d3.axisBottom()
    .scale(xScale);

    svg2.append("g").attr("transform", "translate(0, 100)")
    .call(x_axis);


    


}


function makePopulationData(chromosomeSorted, num){
    var chromDrawarray = {};

    for ( let chromosome of Object.keys(chromosomeSorted)){
        var dataChrom = Object.values(chromosomeSorted[chromosome])
        var m = [].concat.apply([], dataChrom);

        
        var upperLimit = String(Math.max.apply(null, m));
        var currentBase = String(Math.min.apply(null, m));



        var rangeCoord = Math.floor((upperLimit - currentBase)/num)



        chromDrawarray[chromosome]  = new Array(rangeCoord+1).fill(0)

        for ( let coords of dataChrom){
            low = Math.floor((coords[0]-currentBase)/num)
            high= Math.floor((coords[1]-currentBase)/num)

            if (high>chromDrawarray[chromosome].length){

            }
            for (let i = low; i <= high; i++){

                chromDrawarray[chromosome][i]+=1
            }
        }
    }
    return chromDrawarray;

}


function getOrthologs(geneSearched){

    var allOrthologs = new Set();
    allOrthologs.add(geneSearched)
    for ( let gene of orthologs){
        if (gene.source.toLowerCase() == geneSearched.toLowerCase()  || gene.target.toLowerCase() == geneSearched.toLowerCase()){

            allOrthologs.add(gene.source)
            allOrthologs.add(gene.target)
        }
    }
    var data =  Array.from(allOrthologs);
    bringFirst(data, geneSearched)
    return data;
}

function makeMarkers(orthologsArray, chromosomeSorted, geneSearched){

    var pathsToDraw = [];


    removeAll(".markergene")
    for ( let entry of orthologsArray){

        var gene = findGene(entry);


        var svg = document.getElementById(gene.chromosome)

        var chromeStart = getOffset(svg).left-10;
        var svgWidth = svg.getBoundingClientRect().width;


        var markerScale = d3.scaleLinear()
        .domain([0,getChromosomeLength(gene.chromosome,chromosomeSorted)])
    
        .range([0, svgWidth]);

        var markerPos = chromeStart + markerScale(gene.start);

        if(gene.chromosome.slice(0,1)=='D'){

        d3.select(`.genomeView`).append('div').attr("class", "markergene").attr("id",`marker${entry}`).attr("style",`position: absolute;width: 0px;height: 0px;top: 50px;cursor: pointer;opacity: 0.75;left: ${markerPos}; border-bottom: 16px solid orangered; border-left: 10px solid transparent; border-right: 10px solid transparent;`)
        }else{

        d3.select(`.genomeView2`).append('div').attr("class", "markergene").attr("id",`marker${entry}`).attr("style",`position: absolute;width: 0px;height: 0px;top: 50px;cursor: pointer;opacity: 0.75;left: ${markerPos}; border-bottom: 16px solid orangered; border-left: 10px solid transparent; border-right: 10px solid transparent;`)
        
        }
    }
    // https://stackoverflow.com/a/14988898

    if (orthologsArray.length > 1){

    var marker = document.getElementById(`marker${geneSearched}`);
    var sourceGene = {x: getOffset(marker).left+10, y:getOffset(marker).top+8}
    for ( let entry of orthologsArray.slice(1)){
        var marker = document.getElementById(`marker${entry}`);
        pathsToDraw.push({source: sourceGene, target: {x: getOffset(marker).left+10, y:getOffset(marker).top+8 }})

    }

    makePath(pathsToDraw);
}
}


function getOffset(el) {
    const rect = el.getBoundingClientRect();
    return {
      left: rect.left + window.scrollX,
      top: rect.top + window.scrollY
    };
  }

function getChromosomeLength(chromosome, chromosomeSorted){
    var dataChrom = Object.values(chromosomeSorted[chromosome])
    var m = [].concat.apply([], dataChrom);

    
    var upperLimit = String(Math.max.apply(null, m));
    var currentBase = String(Math.min.apply(null, m));



    return(Math.floor(upperLimit - currentBase));
}


function removeAll(className){
    const els = document.querySelectorAll(className);

    els.forEach(el => {
    el.remove();
});
}

const sortAlphaNum = (a, b) => a.localeCompare(b, 'en', { numeric: true })




function createLinkLinePath(d) {
    let curvature = 0.30;
    // code block sourced from d3-sankey https://github.com/d3/d3-sankey for drawing curved blocks
    var x0 = d.source.x,
        x1 = d.target.x,
        y0 = d.source.y,
        y1 = d.target.y,
        yi = d3.interpolateNumber(y0, y1),
        y2 = yi(curvature),
        y3 = yi(1 - curvature);

    return "M" + x0 + "," + y0 + // svg start point
        "C" + x0 + "," + y2 + // curve point 1
        " " + x1 + "," + y3 + // curve point 2
        " " + x1 + "," + y1; // end point
}





function bringFirst(data,toBring){
    data = data.filter(item => item !== toBring);
    data.unshift(toBring);
    }


function makePath(pathsToDraw){

    let sourceGene = pathsToDraw[0]["source"]
    d3.select('#wholeScreen').append("div").attr("id", "paths").attr("style",`position: absolute; top: ${sourceGene["y"]}; left: ${sourceGene["x"]};  width: 100%; height:auto`)
    let screenSVG = d3.select('#paths').append("svg").attr("style",`position: relative; top:0; left:0;width: 100%; height:auto`)
    for (let pair of pathsToDraw){
    
    
    // screenSVG.append("path")
    //         .attr('d', createLinkLinePath({source: {x: 0, y: 0}, target: {x: pair["target"]["x"]-sourceGene["x"], y: pair["target"]["y"]-sourceGene["y"] }}));


    screenSVG.append("line")
    .attr("x1",0)
    .attr("x2",pair["target"]["x"]-sourceGene["x"])
    .attr("y1",0)
    .attr("y2",pair["target"]["y"]-sourceGene["y"])
    .attr("style", "stroke:rgb(255,0,0);stroke-width:2")

}

}