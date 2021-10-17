// sequence metadata class
class SequenceViz
{
    constructor(x, y, width, height, length, hash)
    {
        this.x = x; // locations
        this.y = y;
        this.width = width; //px
        this.height = height; //px
        this.MAXLEN = length;
        this.hash = hash;
    };
};

function main()
{
    const NOD = '3similarhumans' // Name Of Data
    d3.selectAll('title')
        .text(NOD+'.json Visualizer')
        .classed('VizTitle', true)
        .enter();
    d3.select('.header')
        .append('h1')
        .text(NOD+'.json')

    initializeViz(NOD);
}

async function initializeViz(NOD)
{
    const file = 'http://localhost:8000/JSON-OUTPUT/'+NOD+'.json';
    let idata = await d3.json(file);
    vizmain(idata);
}

function vizmain(data)
{
    // init phase
    const library = Object.values(data);
    // console.log(library)
    
    const keys = ["Visible", "Subsequence Length" , "Subsequence Hash"];
    
    // inits html elements of sequences 
    const vizGroup = d3.select('.sequence_box');
    vizGroup.append('div').append('h3').text('NCBI Accession #:')
    const sequenceMeta = Object.values(library[0]);

    const screenMeta = [];
    const renderedSS = [];
    let total_elements = (sequenceMeta.length*2)-1; 

    let c = 0;
    for(let i = 1; i <= total_elements; i++)
    {
        if (i % 2 === 0) 
        { 
            vizGroup.append('svg').classed('lines', true);
        }
        else 
        { 
            // figure out how to center group and rects properly
            let current_class = new SequenceViz(0, 20, 1200
                , 5, sequenceMeta[c].len, sequenceMeta[c].hash);
            screenMeta.push(current_class);
            vizGroup.append('text').text(current_class.hash)
            vizGroup.append('svg')
                .attr('id', 'seq'+formatID(sequenceMeta[c].hash))
                .classed('sequence', true)
                .append('rect')
                .classed('baseline', true)
                .attr('x', current_class.x)
                .attr('y', current_class.y)
                .attr('width', current_class.width)
                .attr('height', current_class.height);
            c++;
        }
    }

    // builds table
    const seqlib = d3.select('.lib')
        .append('table')
        .attr('width', '100%');
    
    const header = seqlib.append('thead');
    header.selectAll('thead').select('th')
        .data(keys)
        .enter()
        .append('th')
        .text(data => data);
    
    const table_data = []; 

    for(let i=1; i<library.length; i++)
    {
        let current = [];
        current.push(library[i].ss_len);
        current.push(library[i].ss_hash);
        table_data.push(current);
    }

    table_data.sort((f, s) => {
        if( s[0] > f[0] ) {return 1}
        else if(f[0] === s[0]) {return 0}
        else{return -1}
    });

    const rows = seqlib.append('tbody');
    rows.selectAll('tr')
        .data(table_data)
        .enter()
        .append('tr')
        .attr('id', (d) => 'row'+d[1]);

    rows.selectAll('tr')
        .append('td')
        .append('input')
        .classed('visbutton', true)
        .attr('type', 'checkbox')
        .attr('id', (d) => d[1])
        .on('change', (d, i) => {updateViewer(i, library, screenMeta, renderedSS, tip);});
        
    const cell = rows.selectAll('tr');
    cell.selectAll('tr')
        .data(data => data)
        .enter()
        .append('td')
        .text(data => data);

    let tip = d3.select('body')
        .append('div')
        .attr('class', 'tooltip')
        .style('opacity', 0);
}

function updateViewer(data, library, screenMeta, renderedSS, tip)
{
    // loads in data and checks for inclusion
    let current_ss = library.filter(d => d.ss_hash === data[1])[0];

    let removed = false;
    let removed_name = '';
    if(renderedSS.filter(d => d.ss_hash == current_ss.ss_hash).length >= 1)
    {   
        let del = renderedSS.findIndex(d => d.ss_hash == current_ss.ss_hash);
        renderedSS.splice(del, 1);
        removed = true;
        removed_name = current_ss.ss_hash;
    }
    else
    {
        renderedSS.push(current_ss);
    }
    let converted_ss = convertSS(renderedSS);

    const current_screen = d3.selectAll('.sequence_box');
    for(e in screenMeta)
    {
        let select_converted_ss = converted_ss.filter( d => {
            if(d[0] == screenMeta[e].hash) {return true}
            else { return false };
        });
        current_screen.selectAll('#seq'+formatID(screenMeta[e].hash))
            .selectAll('.subseq')
            .data(select_converted_ss, d => d)
            .join('rect')
            .classed('subseq', true)
            .attr('id', d => 'rect'+d[1])
            .attr('x', d => {
                return indexToPixel(d[3], screenMeta[e].MAXLEN, screenMeta[e].width);
            })
            .attr('y', screenMeta[e].y-(2*screenMeta[e].height))
            .attr('height', 15)
            .attr('width', d => {
                let ret = indexToPixel(d[2], screenMeta[e].MAXLEN, screenMeta[e].width);
                if(ret === 0)
                {
                    ret = 2;
                };
                return ret
            })
            .attr('fill', d => d[4])
            .merge(current_screen)
            .on('mouseover', (d, i) => {
                // object metadata should contain subsequnece hash as well
                tip.style('opacity', 1)
                    .html(`Subsequence Hash: ${i[1]} <br> Subsequence Length: ${i[2]} bp <br>
                     Start Index: ${i[3]} <br> End Index ${i[3] + i[2]}`);
                d3.selectAll('#row'+i[1])
                    .transition()
                    .duration(200)
                    .style('background-color', i[4]); 
                d3.selectAll('#rect'+i[1])
                    .transition()
                    .duration(50)
                    .style('opacity', 1);
            })
            .on('mousemove', (d) => {
                tip.style('top', d['y']-80+'px')
                    .style('left', d['x']+'px');
            })
            .on('mouseout', (d, i) => {
                tip.style('opacity', 0);
                d3.selectAll('#row'+i[1])
                    .transition()
                    .duration(200)
                    .style('background-color', 'white');
                d3.selectAll('#rect'+i[1])
                    .transition()
                    .duration(50)
                    .style('opacity', .7);
            })
    }
};

function convertSS(ss_list)
{
    let ret = [];
    let localmax = Math.max.apply(null, ss_list.map(d => d.ss_len)); // maybe change to global max later?
    
    for(e in ss_list)
    {
        let ss_len = ss_list[e].ss_len;
        let ss_hash = ss_list[e].ss_hash;
        let current = Object.entries(ss_list[e]).filter(d => {
            if(d[0] === "ss_len" ||  d[0] === "ss_hash") { return false }
            else { return true }
        });
        let color = assignColor(ss_len, localmax);

        // probably a cleaner js solution to this...
        for(i in current)
        {
            for(j in current[i][1])
            {
                // console.log(current[i][1][j])
                // tuple-like object: [main sequence hash, subsequence hash, length of ss in bp, index location in main ss, color]
                ret.push([current[i][0], ss_hash, ss_len, current[i][1][j], color]);
            }
        }
    };
    return ret;
}

function assignColor(ss_len, localmax)
{   
    const hexcolors = ['#ffba08', '#faa307', '#f48c06', '#e85d04', '#dc2f02', '#d00000',];
    let color = Math.floor((ss_len/localmax)*5);
    return hexcolors[color];
}

function formatID(ID)
{
    let formatted = ID.replace(' ', '').replace('.','');
    return formatted;
}

function indexToPixel(index, seqlen, width)
{
    let pixel = Math.floor((index/seqlen)*width);
    if(pixel < 5) { pixel = 5 };
    return pixel;
}

main()