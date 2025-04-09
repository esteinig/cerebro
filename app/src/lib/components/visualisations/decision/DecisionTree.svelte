<script>
    import { onMount } from "svelte";
    import * as d3 from "d3";
  
    // Sample decision tree data structure
    let treeData = {
      question: "Is it raining?",
      answer: "",
      data: { info: "Start of decision tree" },
      children: [
        {
          question: "Do you have an umbrella?",
          answer: "Yes",
          data: { info: "User has an umbrella" },
          children: [
            {
              question: "Will you go outside?",
              answer: "Yes",
              data: { info: "User goes out" },
              children: []
            }
          ]
        },
        {
          question: "Should you stay inside?",
          answer: "No",
          data: { info: "User chooses to stay inside" },
          children: []
        }
      ]
    };
  
    // When the component mounts, create the tree visualization.
    onMount(() => {
      const margin = { top: 20, right: 90, bottom: 30, left: 90 },
        width = 800 - margin.left - margin.right,
        height = 600 - margin.top - margin.bottom;
  
      const svg = d3
        .select("#tree-container")
        .append("svg")
        .attr("width", width + margin.right + margin.left)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
  
      let i = 0,
        duration = 750;
  
      // Create a hierarchy from the data
      let root = d3.hierarchy(treeData, (d) => d.children);
      root.x0 = height / 2;
      root.y0 = 0;
  
      // Optionally collapse nodes after the 1st level for a cleaner initial view
      root.children.forEach(collapse);
  
      function collapse(d) {
        if (d.children) {
          d._children = d.children;
          d._children.forEach(collapse);
          d.children = null;
        }
      }
  
      const treemap = d3.tree().size([height, width]);
  
      update(root);
  
      function update(source) {
        // Compute the new layout.
        const treeDataLayout = treemap(root);
        const nodes = treeDataLayout.descendants(),
          links = treeDataLayout.descendants().slice(1);
  
        // Adjust for fixed-depth
        nodes.forEach((d) => {
          d.y = d.depth * 180;
        });
  
        // **************** Nodes Section ****************
        const node = svg.selectAll("g.node").data(nodes, (d) => d.id || (d.id = ++i));
  
        // Enter any new nodes at the parent's previous position.
        const nodeEnter = node
          .enter()
          .append("g")
          .attr("class", "node")
          .attr("transform", `translate(${source.y0},${source.x0})`)
          .on("click", (event, d) => {
            if (d.children) {
              d._children = d.children;
              d.children = null;
            } else {
              d.children = d._children;
              d._children = null;
            }
            update(d);
          });
  
        // Add circles for nodes.
        nodeEnter
          .append("circle")
          .attr("class", "node-circle")
          .attr("r", 1e-6)
          .style("fill", (d) => (d._children ? "lightsteelblue" : "#fff"));
  
        // Add text labels for nodes.
        nodeEnter
          .append("text")
          .attr("dy", ".35em")
          .attr("x", (d) => (d._children ? -13 : 13))
          .attr("text-anchor", (d) => (d._children ? "end" : "start"))
          .text((d) => d.data.question + (d.data.answer ? ` (${d.data.answer})` : ""))
          .attr("class", "node-text");
  
        // UPDATE nodes.
        const nodeUpdate = nodeEnter.merge(node);
        nodeUpdate
          .transition()
          .duration(duration)
          .attr("transform", (d) => `translate(${d.y},${d.x})`);
  
        nodeUpdate
          .select("circle.node-circle")
          .attr("r", 10)
          .style("fill", (d) => (d._children ? "#555" : "#999"))
          .attr("cursor", "pointer");
  
        // Remove exiting nodes.
        const nodeExit = node
          .exit()
          .transition()
          .duration(duration)
          .attr("transform", (d) => `translate(${source.y},${source.x})`)
          .remove();
  
        nodeExit.select("circle").attr("r", 1e-6);
        nodeExit.select("text").style("fill-opacity", 1e-6);
  
        // **************** Links Section ****************
        const link = svg.selectAll("path.link").data(links, (d) => d.id);
  
        // Enter any new links at the parent's previous position.
        const linkEnter = link
          .enter()
          .insert("path", "g")
          .attr("class", "link")
          .attr("d", (d) => {
            let o = { x: source.x0, y: source.y0 };
            return diagonal(o, o);
          });
  
        // UPDATE links.
        const linkUpdate = linkEnter.merge(link);
        linkUpdate
          .transition()
          .duration(duration)
          .attr("d", (d) => diagonal(d, d.parent));
  
        // Remove any exiting links.
        const linkExit = link
          .exit()
          .transition()
          .duration(duration)
          .attr("d", (d) => {
            let o = { x: source.x, y: source.y };
            return diagonal(o, o);
          })
          .remove();
  
        // Save the old positions for transition.
        nodes.forEach((d) => {
          d.x0 = d.x;
          d.y0 = d.y;
        });
      }
  
      // Creates a curved path from parent to child nodes.
      function diagonal(s, d) {
        return `M ${s.y} ${s.x}
                C ${(s.y + d.y) / 2} ${s.x},
                  ${(s.y + d.y) / 2} ${d.x},
                  ${d.y} ${d.x}`;
      }
    });
  </script>
  
  <style>
    /* Container styling for a vibrant and modern look */
    #tree-container {
      width: 100%;
      height: 600px;
      background: linear-gradient(135deg, #e0f7fa, #80deea);
      border: 1px solid #ccc;
      border-radius: 10px;
      box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
      padding: 10px;
      margin: auto;
    }
  
    /* Node circle style */
    .node-circle {
      stroke: #3182bd;
      stroke-width: 2px;
      transition: stroke-width 0.3s, fill 0.3s;
    }
    
    /* Node text style */
    .node-text {
      font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
      font-size: 12px;
      fill: #333;
      pointer-events: none;
    }
  
    /* Link style */
    .link {
      fill: none;
      stroke: #9ecae1;
      stroke-width: 2px;
    }
  
    /* Interaction hover effects */
    .node:hover circle.node-circle {
      stroke-width: 3px;
      fill: #ffcc00;
    }
  </style>
  
  <div id="tree-container"></div>
  