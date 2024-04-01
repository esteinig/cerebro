
use mongodb::bson::{doc, Document};
use cerebro_model::api::cerebro::model::{SampleId, WorkflowId, CerebroId, RunId, Tag};

/*
============================
MONGDB AGGREGATION PIPELINES
============================
*/


pub fn get_matched_uuid_cerebro_pipeline(id: &CerebroId, taxa: &Option<bool>) -> Vec<Document> {
    
    let pipeline1= vec![
        doc! {
            "$match": {
                "id": &id
            }  
        },
        doc! {
            "$set": {
                "_id": "$$REMOVE",
            }
        },
    ];

    let pipeline2 = vec![
        doc! {
            "$match": {
                "id": &id
            },
        },
        doc! {
            "$set": {
                "_id": "$$REMOVE",
                "taxa": "$$REMOVE",
            }
        },
    ];

    match taxa {
        Some(value) => {
            match value {
                true => pipeline1,
                false => pipeline2,
            }
        },
        None => pipeline1
    }
    
}


pub fn get_matched_id_taxa_cerebro_pipeline(id: &Vec<CerebroId>) -> Vec<Document> {

    vec![
        doc! {
            "$match": {
                "id": { "$in": &id}
            },
        },
        doc! {
            // replace the object root with the nested run doc
            "$replaceRoot": {
                "newRoot": "$taxa"
            }
        },
    ]
}

pub fn get_matched_taxa_summary_pipeline(sample_ids: &Vec<String>, run_ids: &Vec<String>, workflow_ids: &Vec<String>, workflow_names: &Vec<String>) -> Vec<Document> {

    let run_match_stage = match run_ids.is_empty() {
        false => Some(
            doc! {
                "$match": {
                    "run.id": { "$in": &run_ids}
                },
            },
        ),
        true => None
    };

    let sample_match_stage = match sample_ids.is_empty() {
        false => Some(
            doc! {
                "$match": {
                    "sample.id": { "$in": &sample_ids}
                },
            },
        ),
        true => None
    };

    let workflow_match_stage = match workflow_ids.is_empty() {
        false => Some(
            doc! {
                "$match": {
                    "workflow.id": { "$in": &workflow_ids}
                },
            },
        ),
        true => None
    };

    let workflow_name_match_stage = match workflow_names.is_empty() {
        false => Some(
            doc! {
                "$match": {
                    "workflow.name": { "$in": &workflow_names}
                },
            },
        ),
        true => None
    };

    vec![
        workflow_match_stage,
        workflow_name_match_stage,
        run_match_stage,
        sample_match_stage,
        Some(doc! {
            
            "$project": {
                "_id": 0,
                "cerebro_id": "$id",
                "workflow_id": "$workflow.id",
                "workflow_name": "$workflow.name",
                "run_id": "$run.id",
                "run_date": "$run.date",
                "sample_id": "$sample.id",
                "sample_tags": "$sample.tags",
                "sample_type": "$sample.sample_type",
                "sample_group": "$sample.sample_group",
                "taxa": "$taxa"
            }
        }),
    ].into_iter().flatten().collect()  // removes the optional matching stages
}

// pub fn get_matched_id_sample_comments_pipeline(id: &Vec<CerebroId>) -> Vec<Document> {


//     vec![
//         doc! {
//             "$match": {
//                 "id": { "$in": &id}
//             },
//         },
//         doc! {
//             "$group": {
//                 "_id": "$sample.id",
//                 "comments": {
//                     "$addToSet": {
//                         "$unwind": "$sample.comments"
//                     }
//                 }
//             }
//         }
//     ]

// }

pub fn get_paginated_sample_overview_pipeline(page: &i64, limit: &i64, exclude_tags: Vec<&str>, id: &Option<SampleId>, run: &Option<String>, workflow: &Option<String>, group: &Option<String>,) -> Vec<Document> {
    
    let base_match = Some(doc! {
        "$match": {
            "sample.tags": {
                "$nin": exclude_tags
            },
        }
    });

    let id_match: Option<Document> = match id {
        Some(id_query) => Some(doc! {
            "$match": {
                "sample.id": id_query
            }
        }),
        None => None
    };
    let group_match: Option<Document> = match group {
        Some(group_query) => Some(doc! {
            "$match": {
                "sample.group": {
                    "$in": group_query.split(",").map(|x| x.trim()).collect::<Vec<&str>>()
                }
            }
        }),
        None => None
    };

    let run_match: Option<Document> = match run {
        Some(run_query) => Some(doc! {
            "$match": {
                "run.id": {
                    "$in": run_query.split(",").map(|x| x.trim()).collect::<Vec<&str>>()
                }
            }
        }),
        None => None
    };
    let workflow_match = match workflow {
        Some(workflow_query) => Some(doc! {
            "$match": {
                "workflow.id": {
                    "$in": workflow_query.split(",").map(|x| x.trim()).collect::<Vec<&str>>()
                }
            }
        }),
        None => None
    };

    vec![
        base_match,
        id_match,
        group_match,
        run_match,
        workflow_match,

        // Depending on the priority taxa selection (sample.priority) this might get large for memory
        Some(doc! {
            "$project": {
                "run": 1,
                "sample": 1,
                "workflow": 1,
            }
        }),
        Some(doc! {
            // group by sample.id and get overview summaries
            "$group": {
                "_id": "$sample.id",
                "latest_run": {
                    "$max": "$run.date"
                },
                "latest_workflow": {
                    "$max": "$workflow.completed"
                },
                "workflows": {
                    "$addToSet": "$workflow"
                },
                "samples": {
                    "$addToSet": "$sample"
                },
                "runs": {
                    "$addToSet": "$run"
                },
                "description": {
                    "$addToSet": "$sample.description"
                },
                "tags": {
                    "$addToSet": "$sample.tags"
                },
                "groups": { 
                    "$addToSet": "$sample.sample_group"
                },
                "types": { 
                    "$addToSet": "$sample.sample_type"
                },
                "reports": {
                    "$addToSet": "$sample.reports"
                },
                "priority": { 
                    "$addToSet": "$sample.priority" // this should have only the common priority taxa for this sample due to adding to set - it does however not discriminate based on workflow 
                    
                },
            }
        }),
        Some(doc! {
            "$set": {  
                "id": "$_id",
                "_id": "$$REMOVE",
                // Flatten nested fields
                "priority": {
                    "$reduce": {
                        "input": "$priority",
                        "initialValue": [],
                        "in": { "$concatArrays" : ["$$value", "$$this"] }
                    }
                },
                "reports": {
                    "$reduce": {
                        "input": "$reports",
                        "initialValue": [],
                        "in": { "$concatArrays" : ["$$value", "$$this"] }
                    }
                }
            }
        }),
        Some(doc! {
            // sort by sample id
            "$sort": {
                "id": 1
            }
        }),
        Some(doc! {
            // get the page
            "$skip": page*limit
        }),
        Some(doc! {
            // get limited
            "$limit": limit
        }),
        
    ].into_iter().flatten().collect()  // removes the optional matching stages

}

pub fn get_matched_sample_overview_pipeline(id: &str, workflow: &Option<bool>) -> Vec<Document> {
    
    let extract_stage = match workflow {
        Some(true) => {
            doc! {
                "$group": {
                    "_id": "$sample.id",
                    "workflows": { 
                        "$addToSet": "$workflow" 
                    }
                }
            }
        },
        _ => {
            doc! {
                "$group": {
                    "_id": "$sample.id",
                    "data": { 
                        "$addToSet": "$$ROOT"  
                    }
                }
            }
        }
    };

    vec![
        doc! {
            "$match": {
                "sample.id": &id
            }
        },
        doc! {
            "$project": {
                "run": 1,
                "sample": 1,
                "workflow": 1,
            }
        },
        extract_stage,
        doc! {
            "$set": {
                "id": "$_id",
                "_id": "$$REMOVE",
            }
        }
    ]

}


pub fn get_matched_sample_cerebro_notaxa_pipeline(id: &str, workflow: &Option<WorkflowId>) -> Vec<Document> {
    
    let match_stage = match workflow {
        Some(wf_id) => {
            doc! {
                "$match": {
                    "sample.id": &id,
                    "workflow.id": &wf_id
                }
            }
        },
        _ => {
            doc! {
                "$match": {
                    "sample.id": &id
                }
            }
        }
    };

    vec![
        match_stage,
        doc! {
            "$set": {
                "_id": "$$REMOVE",
                "taxa": "$$REMOVE"
            }
        }
    ]
}

pub fn get_matched_samples_cerebro_notaxa_pipeline(sample_ids: &Vec<String>, cerebro_ids: &Vec<String>, workflow: &Option<WorkflowId>) -> Vec<Document> {
    
    let match_stage = match workflow {
        Some(wf_id) => {
            match (sample_ids.is_empty(), cerebro_ids.is_empty()) {
                (false, true) => {
                    doc! {
                        "$match": {
                            "sample.id": { "$in": &sample_ids },
                            "workflow.id": &wf_id
                        }
                    }
                },
                (true, false) => {
                    doc! {
                        "$match": {
                            "id": { "$in": &cerebro_ids },
                            "workflow.id": &wf_id
                        }
                    }
                },
                // If neither, return empty aggregate pipeline
                (true, true) => {
                    return vec![]
                },
                // If both use more specific model identifiers
                (false, false) => {
                    doc! {
                        "$match": {
                            "id": { "$in": &sample_ids },
                            "workflow.id": &wf_id
                        }
                    }
                }
            }
        },
        _ => {
            match (sample_ids.is_empty(), cerebro_ids.is_empty()) {
                (false, true) => {
                    doc! {
                        "$match": {
                            "sample.id": { "$in": &sample_ids }
                        }
                    }
                },
                (true, false) => {
                    doc! {
                        "$match": {
                            "id": { "$in": &cerebro_ids }
                        }
                    }
                },
                // If neither, return empty aggregate pipeline
                (true, true) => {
                    return vec![]
                },
                // If both use more specific model identifiers
                (false, false) => {
                    doc! {
                        "$match": {
                            "id": { "$in": &sample_ids }
                        }
                    }
                }
            }
        }
    };

    vec![
        match_stage,
        doc! {
            "$project": {
                "_id": 0,
                "id": "$id",
                "quality": "$quality",
                "run": "$run",
                "sample": "$sample",
                "workflow": "$workflow"
            }
        }
    ]
}

pub fn get_matched_sample_ids_cerebro_pipeline(ids: &Vec<String>, workflow: &Option<WorkflowId>) -> Vec<Document> {
    
    let match_stage = match workflow {
        Some(wf_id) => {
            doc! {
                "$match": {
                    "sample.id": { "$in": &ids },
                    "workflow.id": &wf_id
                }
            }
        },
        _ => {
            doc! {
                "$match": {
                    "sample.id": { "$in": &ids }
                }
            }
        }
    };

    vec![
        match_stage,
        doc! {
            "$project": {
                "_id": 0,
                "sample_id": "$sample.id"
            }
        }
    ]
}


pub fn get_matched_workflow_cerebro_notaxa_pipeline(workflow: &str, runs: &Option<String>, tags: &Option<String>,  return_uuid: &Option<bool>) -> Vec<Document> {
    
    let default = vec![
        doc! {
            "$match": {
                "workflow.id": &workflow
            }
        },
        doc! {
            "$set": {
                "_id": "$$REMOVE"
            }
        },
    ];

    let mut pipeline = match tags {
        Some(tag_str) => {
            let include_tags = tag_str.split(",").map(|x| x.trim().to_string()).collect::<Vec<Tag>>();
                vec![
                    doc! {
                        "$match": {
                            "sample.tags": {
                                "$in": include_tags
                            }
                        }
                    }
                ]
        },
        None => default
    };

    if let Some(run_str) = runs {
        let include_runs = run_str.split(",").map(|x| x.trim().to_string()).collect::<Vec<RunId>>();
            pipeline.push(
                doc! {
                    "$match": {
                        "run.id" : {
                            "$in": include_runs
                        }
                    }
                }
            )
    };

    if let Some(return_uuid) = return_uuid {
        if *return_uuid {
            pipeline.push(doc! {
                "$project": { "id": 1, "_id": 0 }
            })
        }
    }

    pipeline
}