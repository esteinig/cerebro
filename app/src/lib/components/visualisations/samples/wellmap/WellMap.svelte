<script lang="ts">

    import { onMount } from 'svelte';
    import * as d3 from 'd3';
  
    // Data should be an array of objects: { id: string, value: number, score: number }
    export let data: any[] = [
        {"id": "V01", "value": 5.0, "score": 3},
        {"id": "V02", "value": 1.0, "score": 3},
        {"id": "V03", "value": 100.0, "score": 3},
        {"id": "V04", "value": 10.0, "score": 3},
        {"id": "V05", "value": 15.0, "score": 3},
        {"id": "V06", "value": 5.0, "score": 3},
    ];
    export let numCol = 12;
    export let numRow = 8;

    // Default color scale for filling the well circles (customisable)
    export let colorScale = d3.scaleSequential(d3.interpolateViridis);
    
    // Default ordinal scale for the outer score shells (customisable)
    export let shellScale = d3.scaleOrdinal(d3.schemeCategory10);
  
    // Adjustable cell size for each grid position
    export let cellSize = 100;
  
    let svg: any;

    // Define margins and based on the grid and cell size
    const margin = { top: 20, right: 20, bottom: 20, left: 20 };

    const width = numCol * cellSize + margin.left + margin.right;
    const height = numRow * cellSize + margin.top + margin.bottom;
    
    // Map the input data into grid cells based on row-major order
    const cells: any[] = [];

    for (let r = 0; r < numRow; r++) {
        for (let c = 0; c < numCol; c++) {
            const idx = r * numCol + c;
            if (idx < data.length) {
                cells.push({ ...data[idx], row: r, col: c });
            }
        }
    }

    // Define the radius for the main well circle and the shell thickness for rings.
    const circleRadius = cellSize * 0.3;
    const shellThickness = 3;
    const shellGap = 2;

    const range = (n: number) => Array.from({ length: n }, (_, i) => i + 1);

  </script>
  
  <svg bind:this={svg} width={width} height={height}>
    <g transform="translate({margin.left},{margin.top})">
        {#each cells as cell}
            <g class="cell" transform="translate({cell.col * cellSize + cellSize / 2}, {cell.row * cellSize + cellSize / 2})">
                <circle r={circleRadius} fill={colorScale(cell.value)}>
                    
                </circle>
                {#each range(cell.score) as i}
                    <circle r={circleRadius + shellGap + shellThickness/2 + (i - 1) * (shellGap + shellThickness)} fill="none" stroke={shellScale(i)} stroke-width={shellThickness}></circle>
                {/each}
            </g>
        {/each}
    </g>
  </svg>
  