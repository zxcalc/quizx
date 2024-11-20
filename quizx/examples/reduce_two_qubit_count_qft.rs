use quizx::circuit::*;
use quizx::vec_graph::*;
use quizx::approximate::*;

fn main() {
    let f = "circuits/small/qft_8.qasm";
    println!("{}", f);
    let err_array = vec![0.0,0.01111111,0.02222222,0.03333333,0.04444444,0.05555556,0.06666667,0.07777778,0.08888889,0.1];
    let mut anneal_gate_mean_array = vec![0.0; 10];
    let mut anneal_gate_min_array = vec![0.0; 10];
    let mut greedy_gate_array = vec![0; 10];

    let c = Circuit::from_file(f)
            .unwrap_or_else(|_| panic!("circuit failed to parse: {}", f))
            .to_basic_gates();
    let g: Graph = c.to_graph();
    let inital_count = two_qubit_count(&g);

    // Iterate through each error budget
    for (i,&err_budget) in err_array.iter().enumerate() {
        let mut iter_for_anneal = vec![0;4];
        // Iterate a few times to collect mean and minimum
        for j in 0..4 {
            // Greedy approximation
            let c = Circuit::from_file(f)
                .unwrap_or_else(|_| panic!("circuit failed to parse: {}", f))
                .to_basic_gates();
            let g_out = basic_greedy(c, err_budget);
            greedy_gate_array[i] = two_qubit_count(&g_out);
            
            // Annealed approximation
            let c = Circuit::from_file(f)
                .unwrap_or_else(|_| panic!("circuit failed to parse: {}", f))
                .to_basic_gates();
            let g: Graph = c.to_graph();
            let g_out: Graph = two_qubit_basic_anneal(g, err_budget, 100000,0.9, true);
            iter_for_anneal[j] = two_qubit_count(&g_out);
        }
        let mean = iter_for_anneal.iter().sum::<usize>() as f32 / iter_for_anneal.len() as f32;
        anneal_gate_mean_array[i] = mean;
        let min_value = iter_for_anneal.iter().min().copied().unwrap_or(0) as f64;
        anneal_gate_min_array[i] = min_value;
    }
    println!("Error array: {}",err_array.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    println!("Initial count: {}",inital_count.to_string().as_str());
    println!("Annnealed gate count (mean): {}",anneal_gate_mean_array.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    println!("Annnealed gate count (min): {}",anneal_gate_min_array.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
    println!("Greedy gate count: {}",greedy_gate_array.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
}