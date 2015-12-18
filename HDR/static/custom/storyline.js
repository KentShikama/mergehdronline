$(document).ready(function () {
    var fiction_id = $('#fiction_id').text();

    var min_size = 20; var max_size = 200; var proportion = 10;

    function calculateColorBasedOnSentiment(sentiment){
        var adjusted_sentiment = Math.floor((parseFloat(sentiment) + 0.99) * 5);
        var colors = ["#D4253F", "#D24525", "#D07D26", "#CEB426", "#AECC26", "#76CA27", "#40C827", "#28C644", "#28C477", "#29C2AA"];
        return colors[adjusted_sentiment];
    }

    function episodeInformation(id, summary) {
        return summary +
            "<div class='margin20'></div>" +
            "<a href='/episode/" + id + "'><div class='btn btn-default'>Expand in episode view</div></a>" +
            "<div class='margin10'></div>" +
            "<a href='/episode/new/" + fiction_id + "/" + id + "'><div class='btn btn-default'>Evolve this episode</div></a>"+
            "<div class='margin10'></div>" +
            "<a href='/storyline/" + id + "/" +  "'><div class='btn btn-default'>View Storyline from Beginning</div></a>";
    }

    var cy = window.cy = cytoscape({
                                       container: document.getElementById('cy'),

                                       boxSelectionEnabled: false,
                                       autounselectify: true,

                                       layout: {
                                           name: 'dagre'
                                       },

                                       style: [
                                           {
                                               selector: 'node',
                                               style: {
                                                   'content': 'data(label)',
                                                   'text-valign': 'center',
                                                   'text-halign': 'right',
                                                   'background-color': 'data(color)',
                                                   'color': '#333',
                                                   'text-wrap': 'wrap',
                                                   'text-max-width': 160,
                                                   'shape': 'circle',
                                                   'width': 'data(popularity)',
                                                   'height': 'data(popularity)'
                                               }
                                           },

                                           {
                                               selector: 'edge',
                                               style: {
                                                   'width': 4,
                                                   'opacity': 0.4,
                                                   'target-arrow-shape': 'triangle',
                                                   'line-color': '#666',
                                                   'target-arrow-color': '#666'
                                               }
                                           }
                                       ]
                                   });

    function handle_node_hover() {
        cy.on('mouseover', 'node', function (event) {
            var node = event.cyTarget;
            var summary = node.attr('summary');
            node.qtip({
                          content: summary,
                          position: {
                              at: 'bottom center'
                          },
                          show: {
                              event: "mouseover",
                              solo: true,
                              ready: true
                          },
                          hide: {
                              event: 'mouseout'
                          },
                          style: {
                              classes: 'qtip-bootstrap',
                              tip: {
                                  width: 16,
                                  height: 8
                              }
                          }
                      }, event);
        });
    }

    function handle_nodes_creation(content) {
        var edge_list = [];
        jsonQ.each(content, function (key, value) {
            var node = {group: "nodes", data: {}};
            var id = value.id;
            node.data["id"] = value.id;
            node.data["label"] = value.title;
            node.data["color"] = calculateColorBasedOnSentiment(value.sentiment);
            node.data["summary"] = episodeInformation(id, value.summary);
            node.data["popularity"] = Math.min(min_size+proportion*Math.log2(value.popularity+2),max_size);
            cy.add(node);
            jsonQ.each(value.previous_ids_without_parent, function (key, value) {
                var edge = {group: "edges", data: {}};
                edge.data["source"] = value;
                edge.data["target"] = id;
                edge_list.push(edge);
            });
            jsonQ.each(value.next_ids_without_parent, function (key, value) {
                var edge = {group: "edges", data: {}};
                edge.data["source"] = id;
                edge.data["target"] = value;
                edge_list.push(edge);
            });
        });
        jsonQ.each(edge_list, function (key, value) {
            cy.add(value);
        });
    }

    $.ajax({url: "/api/fiction/" + fiction_id}).then(function (content) {
        handle_nodes_creation(content);
        cy.layout({name: "dagre"});
        handle_node_hover();
    }, function (xhr, status, error) {
        console.log(error);
    });
});