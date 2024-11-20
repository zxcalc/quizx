use crate::circuit::*;
use crate::vec_graph::*;
use rand::prelude::*;
use std::collections::HashMap;
use crate::simplify::basic_simp;

pub fn basic_greedy(c: Circuit, err_budget: f64) -> Graph {
    let mut g: Graph = c.to_graph();
    let mut count = two_qubit_count(&g);
    println!("Two-qubit count before: {count}");

    let mut total_err: f64 = 0.0;

    let squashable_vertices: Vec<usize> = g.vertices()
        .into_iter()
        .filter(|&i| g.phase(i).to_f64() % 0.5 != 0.0)
        .collect();

    for v in squashable_vertices{
        let phase = g.phase(v).to_f64();
        let phase_rounded = (phase * 2.0).round()/2.0;
        let difference = phase_rounded - phase;
        let norm_err = (((difference/2.0) * std::f64::consts::PI).sin()).abs() * 2.0;
        if (total_err + norm_err) < err_budget{
            let mut g_test = g.clone();
            g_test.set_phase(v,phase_rounded);
            basic_simp(&mut g_test);
            let count_new = two_qubit_count(&g_test);
            if count_new < count{
                count = count_new;
                g.set_phase(v,phase_rounded);
                total_err += norm_err;
            };
        };
    };
    println!("Total error: {}", total_err);
    basic_simp(&mut g);
    let count_after = two_qubit_count(&g);
    println!("Two-qubit count after: {count_after}");
    g
}

pub fn two_qubit_basic_anneal(mut g: Graph, err_budget: f64, max_iter: usize, cooling_rate: f64, quiet:bool) -> Graph {
    // General parameters
    let mut total_err: f64 = 0.0;
    let mut tq = two_qubit_count(&g);
    let mut tq_best = tq;
    let mut rng = thread_rng();
    let mut g_best = g.clone();
    let reason_for_break:&str; 

    let mut change_counter = 0;
    let mut no_change_counter = 0;
    let mut freeze_out_counter = 0;
    let mut iter = 0;
    let mut new_best_reached = false;

    let mut temp = 1.0;
    let epoch = g.num_vertices();
    let mut temp_array = vec![temp];

    // Collect non-Clifford vertices
    let mut vertex_list: Vec<usize> = g
        .vertices()
        .filter(|&v| g.phase(v).to_f64() % 0.5 != 0.0)
        .collect();

    let mut vertex_flip: HashMap<usize, bool> = vertex_list.iter().map(|&v| (v, false)).collect();
    let phase_dict: HashMap<usize, f64> = vertex_list.iter().map(|&v| (v, g.phase(v).to_f64())).collect();
    let rounded_phase_dict: HashMap<usize, f64> = vertex_list
        .iter()
        .map(|&v| (v, (2.0 * g.phase(v).to_f64()).round() / 2.0))
        .collect();
    let mut round_error_dict: HashMap<usize, f64> = HashMap::new();
    for &v in &vertex_list {
        let phase = phase_dict[&v];
        let rounded_phase = rounded_phase_dict[&v];
        let error = (rounded_phase - phase).abs();
        round_error_dict.insert(v, error);
    }

    let freeze_out_max= vertex_list.len() * 3;

    // Remove vertices that would exceed the error budget (rough upper bound)
    vertex_list.retain(|&v| total_err + (round_error_dict[&v] * 1.5) <= err_budget);
    if vertex_list.is_empty() {
        println!("All gates would break error budget, nothing to squash.");
        return g;
    }
    // Create probability list
    let inv_error: Vec<f64> = vertex_list.iter().map(|&v| 1.0 / round_error_dict[&v]).collect();
    let prob_sum: f64 = inv_error.iter().sum();
    let prob_to_be_picked: Vec<f64> = inv_error.iter().map(|&e| e / prob_sum).collect();

    while total_err <= err_budget && iter < max_iter && freeze_out_counter < freeze_out_max && tq != 0 {
        let mut iter_at_temp = 0;

        // Loop at a given temperature till break conditions
        loop {
            // Pick a random vertex based on probability
            let v = *vertex_list
                .choose_weighted(&mut rng, |&i| prob_to_be_picked[vertex_list.iter().position(|&x| x == i).unwrap()])
                .unwrap();

            let norm_err = 2.0 * (round_error_dict[&v] * std::f64::consts::PI / 2.0).sin().abs();
            // If the vertex is not rounded and error would break the budget, skip
            if !vertex_flip[&v] && (total_err + norm_err) > err_budget {
                freeze_out_counter += 1;
            } else {
                // Flip the vertex
                let mut g_test = g.clone();
                if !vertex_flip[&v]{
                    g_test.set_phase(v, rounded_phase_dict[&v]);
                } else {
                    g_test.set_phase(v, phase_dict[&v]);
                }
                basic_simp(&mut g_test);
                let tq_new = two_qubit_count(&g_test);

                // Better result acceptance
                if tq_new < tq {
                    freeze_out_counter = 0;
                    change_counter += 1;
                    // Change the error
                    if !vertex_flip[&v] {
                        g.set_phase(v, rounded_phase_dict[&v]);
                        total_err += norm_err;
                    } else {
                        g.set_phase(v, phase_dict[&v]);
                        total_err -= norm_err;
                    }
                    // Change the new tq count
                    tq = tq_new;
                    // Log that the vertex has been toggled
                    vertex_flip.insert(v, !vertex_flip[&v]);
                    if tq_new < tq_best {
                        tq_best = tq_new;
                        new_best_reached = true;
                        g_best = g.clone();
                    }
                    if tq_best == 0 {
                        break;
                    }
                } else if tq_new > tq && rng.gen::<f64>() < acceptance_probability(tq_new, tq, temp) {
                    freeze_out_counter = 0;
                    change_counter += 1;
                    // Change the error
                    if !vertex_flip[&v] {
                        g.set_phase(v, rounded_phase_dict[&v]);
                        total_err += norm_err;
                    } else {
                        g.set_phase(v, phase_dict[&v]);
                        total_err -= norm_err;
                    }
                    tq = tq_new;
                    vertex_flip.insert(v, !vertex_flip[&v]);
                } else {
                    vertex_flip.insert(v, !vertex_flip[&v]);
                    no_change_counter += 1;
                    freeze_out_counter += 1;
                }
            }
            iter += 1;
            iter_at_temp += 1;
            temp_array.push(temp);

            // Check if new best has been found
            if iter_at_temp % epoch == 0 && !new_best_reached{
                new_best_reached = false;
                break;
            }
            else if iter_at_temp % epoch == 0 {
                new_best_reached = false;
            }

            if freeze_out_counter >= freeze_out_max || (iter >= max_iter) {
                break;
            }
        }
        temp *= cooling_rate;
    }
    if iter >= max_iter {
        reason_for_break = "Max iterations exceeded";
    } else if tq == 0 {
        reason_for_break = "Two-qubit count is zero";
    } else if freeze_out_counter == freeze_out_max {
        reason_for_break = "System froze";
    } else {
        reason_for_break = "Unknown break";
    }

    if !quiet {
        println!("{}", reason_for_break);
        println!("Changes / Not Changes: {change_counter} / {no_change_counter}");
        println!("Total error: {}", total_err);
        println!("Two-qubit count after annealing: {}", tq_best);
    }
    basic_simp(&mut g_best);
    g_best
}

fn acceptance_probability(new_tq: usize, old_tq: usize, temp: f64) -> f64 {
    (-((new_tq - old_tq) as f64) / temp).exp()
}

pub fn two_qubit_count(g: &impl GraphLike) -> usize {
    g.num_edges() - g.num_vertices() + g.inputs().len()
}