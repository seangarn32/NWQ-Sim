#include <vector>
#include <string>
#include <map>
#include <complex>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <map>
#include <exception>
#include <utility>
#include <cmath>
#include "pcg_random.hpp"
#include "../../Eigen/Dense"

//Using dynamic Eigen MatrixXd to streamline calculations
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix2D = Eigen::Matrix2cd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;

typedef std::pair<std::string, double> Term; //A fermionic string <operator, coefficient>
typedef std::vector<Term> CircuitConfig; //In between circuit for translation from fermionic operators to pauli operators

int generateRandomNumber(uint64_t seed) {
    //Create a PCG random number generator initialized with the given seed
    pcg64 rng(seed);
    std::uniform_int_distribution<int> dist(0, 100); // Distribution that maps to the range [0,100]
    return dist(rng); // Generate and return the random number
}


//Generates the t values in the correct order depending if its the left or right operator
void Gen_tLR(Vector& tLR, const std::string& name, const Vector& ts, bool left) { 
    size_t count_R = std::count(name.begin(), name.end(), 'R');
    //Single operator case
    if (count_R == 1) {
        tLR[0] = ts[std::stoi(name.substr(name.find('R') + 1))];
    //Composite case
    } else if (count_R == 2) {
        std::string filtered = std::regex_replace(name, std::regex("[+-]"), "");
        std::vector<int> indices;
        size_t pos = 0;
        while ((pos = filtered.find('R', pos)) != std::string::npos) {
            indices.push_back(ts[std::stoi(filtered.substr(pos + 1, filtered.find_first_not_of("0123456789", pos + 1) - pos - 1))]);
            pos += 2; //Move past 'R' and the digit
        }
        if (indices.size() == 2) {
            //If its the left bases
            if(left)
                tLR = {indices[1], indices[0]};
            //If its the right bases
            else
                tLR = {indices[0], indices[1]};
        } else {
            throw std::runtime_error("Invalid basis name: " + name);
        }
    } else {
        throw std::runtime_error("Invalid basis name: " + name);
    }
}

class PauliOperator {
public:
    std::map<std::string, double> terms;

    void addTerm(const std::string& term, double coefficient) {
        if (terms.find(term) == terms.end()) {
            terms[term] = coefficient;
        } else {
            terms[term] += coefficient;
        }
    }

    void display() const {
        for (const auto& term : terms) {
            std::cout << term.first << " * " << term.second << "\n";
        }
    }
};

//Functions to define individual Pauli matrices
ComplexMatrix2D pauliI() {
    ComplexMatrix2D I;
    I << 1, 0,
         0, 1;
    return I;
}

ComplexMatrix2D pauliX() {
    ComplexMatrix2D X;
    X << 0, 1,
         1, 0;
    return X;
}

ComplexMatrix2D pauliY() {
    ComplexMatrix2D Y;
    Y << 0, -std::complex<double>(0, 1),
         std::complex<double>(0, 1), 0;
    return Y;
}

ComplexMatrix2D pauliZ() {
    ComplexMatrix2D Z;
    Z << 1, 0,
         0, -1;
    return Z;
}

//Maps the Pauli matrices to the PauliOps
std::map<char, ComplexMatrix2D> createPauliMap() {
    std::map<char, ComplexMatrix2D> pauliMap;
    pauliMap['I'] = pauliI();
    pauliMap['X'] = pauliX();
    pauliMap['Y'] = pauliY();
    pauliMap['Z'] = pauliZ();
    return pauliMap;
}

//Converts a single Pauli string to a matrix
ComplexMatrix pauliStringToMatrix(const std::string& pauliString) {
    auto pauliMap = createPauliMap();
    size_t numQubits = pauliString.length();
    ComplexMatrix result = ComplexMatrix::Identity(1, 1);

    for (size_t i = 0; i < numQubits; i++) {
        ComplexMatrix2D currentMatrix = pauliMap[pauliString[i]];
        result = Eigen::kroneckerProduct(result, currentMatrix).eval();
    }

    return result;
}

//Converts a PauliOperator to a full matrix
ComplexMatrix pauliOperatorToMatrix(const PauliOperator& pauliOp, int numQubits) {
    size_t matrixDim = 1 << numQubits; // 2^numQubits
    ComplexMatrix result = ComplexMatrix::Zero(matrixDim, matrixDim);

    for (const auto& term : pauliOp.terms) {
        const std::string& pauliString = term.first;
        double coefficient = term.second;
        result += coefficient * pauliStringToMatrix(pauliString);
    }

    return result;
}

//Maps fermionic creation and annihilation operators to Pauli strings
std::pair<std::string, std::string> fermionicToPauli(int index, int numQubits) {
    std::string creation(numQubits, 'I');
    std::string annihilation(numQubits, 'I');

    //Add the Z strings up to index
    for (int i = 0; i < index; i++) {
        creation[i] = 'Z';
        annihilation[i] = 'Z';
    }

    //Add X and Y Pauli operators at the index position
    creation[index] = 'X';
    annihilation[index] = 'Y';

    return std::make_pair(creation, annihilation);
}

//Inverted Jordan-Wigner Transformation Function
PauliOperator jordanWignerTransform(const CircuitConfig& circuitConfig, int numQubits) {
    PauliOperator pauliOperator;
    for (const auto& fop : circuitConfig) {
        std::string fermionStr = fop.first;
        double coefficient = fop.second;

        for (size_t i = 0; i < fermionStr.length(); i++) {
            if (isdigit(fermionStr[i])) {
                int index = fermionStr[i] - '0';
                if (index < 0 || index >= numQubits) {
                    throw std::out_of_range("Index out of range for the number of qubits");
                }

                bool isCreation = fermionStr[i - 1] == '+';
                std::pair<std::string, std::string> pauliStrings = fermionicToPauli(index, numQubits);

                pauliOperator.addTerm(isCreation ? pauliStrings.first : pauliStrings.second, coefficient);
            }
        }
    }
    return pauliOperator;
}

//Processes a single R into a fermionic string
CircuitConfig baseCircuitCodeSingle(const std::string& singleName, double t,
                                    const std::map<std::string, CircuitConfig>& basisDict, bool left) {
    CircuitConfig totalConfig;
    char signStr = singleName[0]; // '+' or '-'
    double tAct;

    if (signStr == '+') {
        tAct = left ? -t : t;
    } else if (signStr == '-') {
        tAct = left ? t : -t;
    } else {
        throw std::invalid_argument("Unknown sign in singleName");
    }

    std::string baseStr = singleName.substr(1); // 'R<num>', e.g., 'R1'
    auto baseConfigIter = basisDict.find(baseStr);
    if (baseConfigIter == basisDict.end()) {
        throw std::invalid_argument("Unknown basis string: " + baseStr);
    }

    const CircuitConfig& baseConfig = baseConfigIter->second;
    for (const Term& term : baseConfig) {
        totalConfig.emplace_back(term.first, term.second * tAct);
    }

    return totalConfig;
}

//Exponentiates a single fermionic string and return its matrix representation
ComplexMatrix expSingleR(const std::string& singleName, double t,
                                  const std::map<std::string, CircuitConfig>& basisDict,
                                  int numQubits, bool left) {
    CircuitConfig config = baseCircuitCodeSingle(singleName, t, basisDict, left);
    PauliOperator pauliOp = jordanWignerTransform(config, numQubits);
    ComplexMatrix H = pauliOperatorToMatrix(pauliOp, numQubits);
    return (-std::complex<double>(0, 1) * H).exp(); // e^{-i * H}
}

//Parses composite fermionic strings and compute individual exponentials in left to right order
ComplexMatrix expCompositeRight(const std::string& nameComposite, const Vector& ts,
                                         const std::map<std::string, CircuitConfig>& basisDict,
                                         int numQubits, bool left) {
    std::vector<std::string> subNames;
    Vector timeParams = ts;

    //Splits the composite name by "+" or "-" signs while keeping the delimiters
    size_t pos = 0, found;
    while ((found = nameComposite.find_first_of("+-", pos)) != std::string::npos) {
        if (found > pos) {
            subNames.push_back(nameComposite.substr(pos, found - pos));
        }
        subNames.push_back(std::string(1, nameComposite[found]) + nameComposite.substr(found + 1, nameComposite.find_first_of("+-", found + 1) - (found + 1)));
        pos = found + 1;
    }
    subNames.push_back(nameComposite.substr(pos));

    //Adjusts the signs based on the left parameter
    for (size_t i = 0; i < subNames.size(); i++) {
        timeParams[i] *= (subNames[i][0] == '+' ? 1.0 : -1.0) * (left ? -1.0 : 1.0);
    }

    //Computes the exponential of each individual sub-operator sequentially left to right
    ComplexMatrix compositeExp = ComplexMatrix::Identity(1 << numQubits, 1 << numQubits);
    for (size_t i = 0; i < subNames.size(); i++) {
        compositeExp *= expSingleR(subNames[i], timeParams[i], basisDict, numQubits, left);
    }

    return compositeExp;
}

//Parses composite fermionic strings and compute individual exponentials in right to left order
ComplexMatrix expCompositeLeft(const std::string& nameComposite, const Vector& ts,
                                        const std::map<std::string, CircuitConfig>& basisDict,
                                        int numQubits, bool left) {
    std::vector<std::string> subNames;
    Vector timeParams = ts;

    //Splits the composite name by "+" or "-" signs while keeping the delimiters
    size_t pos = 0, found;
    while ((found = nameComposite.find_first_of("+-", pos)) != std::string::npos) {
        if (found > pos) {
            subNames.push_back(nameComposite.substr(pos, found - pos));
        }
        subNames.push_back(std::string(1, nameComposite[found]) + nameComposite.substr(found + 1, nameComposite.find_first_of("+-", found + 1) - (found + 1)));
        pos = found++;
    }
    subNames.push_back(nameComposite.substr(pos));

    //Adjusts the signs based on the left parameter
    for (size_t i = 0; i < subNames.size(); i++) {
        timeParams[i] *= (subNames[i][0] == '+' ? 1.0 : -1.0) * (left ? -1.0 : 1.0);
    }

    //Computes the exponential of each individual sub-operator sequentially right to left
    ComplexMatrix compositeExp = ComplexMatrix::Identity(1 << numQubits, 1 << numQubits);
    for (int i = subNames.size() - 1; i >= 0; --i) {
        compositeExp *= expSingleR(subNames[i], timeParams[i], basisDict, numQubits, left);
    }

    return compositeExp;
}

//Computes matrix elements classically
void ComputeSH(ComplexMatrix HFState, ComplexMatrix HamMat, Matrix& S, Matrix& H, const std::vector<std::string>& bases, const Vector& ts, 
               const std::map<std::string, CircuitConfig>& basisDict, int numOrbitals, int matLen, bool debug = false) {
    //Start with the left side of <R| |R> as the outer loop to fill in the matrices
    //Initialize some strings for the name of the right and left bases
    std::string nar, nal;
    //Initialize some intermediate matrices
    ComplexMatrix HMatMeasurement, HMatMeasurable, HMatExpectation, SMatMeasurement, SMatMeasurable, SMatExpectation;
    ComplexMatrix SObs, HObs;
    //Initialize the final complex matrices, the real part wil be place in S and H
    ComplexMatrix SComp, HComp;
    //The t's that will correspond with the right and left bases respectively
    Vector tr, tl;
    //Left and right pauli operators and matrix representations of the Pauli operators
    ComplexMatrix leftMatrix, rightMatrix;
    //Outer loop & left basis(bases) calculation
    for (int i = 0; i < matLen; i++) {
        nal = bases[i]; //Name left
        tl(-1, -1); //Default to -1 no value
        Gen_tLR(tl, nal, ts, true); //Get the pair of t's in order for the left side
        //Output the tuple
        //std::cout << "Tuple: (" << tl.first << ", " << (tl.second == -1 ? "" : std::to_string(tl.second)) << ")\n";
        //Generate pauli operators using the basis and basis dictionary and perform JW transform
        //Convert the operators into a matrix, each orbital corresponds to a qubit
        leftMatrix = expCompositeLeft(nal, tl, basisDict, numOrbitals, debug);
        
        //Take the right side of <R| |R> as the inner loop to fill in the matrices
        //Inner loop & right basis(bases) calculation
        for (int j = 0; j < matLen; j++) {
            //Name of right basis(bases)
            nar = bases[j];
            //t values that we are interested in for this(these) particular basis(bases)
            tr(-1, -1);
            //Get the t or pair of t's in order for the right side
            Gen_tLR(tr, nar, ts, false);

            HMatMeasurement = HamMat.adjoint(); 
            //Diagonal entries in S and H matrices
            if(nal == nar){ 
                SComp[i,j] = 1;
                //Simplify a bit in R0 basis case
                if(nal == "R0"){
                    HMatMeasurable = HMatMeasurement  * HFState;
                }
                else{
                    rightMatrix = expCompositeRight(nar, tr, basisDict, numOrbitals, debug);
                    HMatMeasurable = HMatMeasurement * rightMatrix * HFState;
                }
                HComp[i,j] = HMatMeasurable.eval();
            }
            //Off diagonal
            else{ 
                //Simplify a bit in R0 basis case
                if(nar == "R0"){
                    SObs = leftMatrix;
                    HObs = leftMatrix * HamMat;
                }
                //Simplify a bit in R0 basis case
                else if (nal == "R0"){
                    rightMatrix = expCompositeRight(nar, tr, basisDict, numOrbitals, debug);
                    SObs = rightMatrix;
                    HObs = HamMat * rightMatrix;
                }
                else{
                    rightMatrix = expCompositeRight(nar, tr, basisDict, numOrbitals, debug);
                    SObs = leftMatrix * rightMatrix;
                    HObs = leftMatrix * HamMat * rightMatrix;
                }
                //S matrix
                SMatMeasurement = SObs.adjoint();
                SMatMeasurable = SMatMeasurement * HFState;
                SComp[i,j] =  SMatMeasurable.eval();
                //H matrix
                HMatMeasurement = HObs.adjoint();
                HMatMeasurable = HMatMeasurement * HFState;
                HComp[i,j] = HMatMeasurable.eval();
            }
        }
    }
    S = SComp.real();
    H = HComp.real();
}

//for each Vp;
//Step 7 & 8
//Step 9 & 10: Evaluate matrix elements on a quantum device
//Matrix evaluatedMatrixElements = EvaluateMatrixElements(matrixElements);
//Function to evaluate matrix elements on a quantum device

void QuantumSH(ComplexMatrix HFState, ComplexMatrix HamMat, Matrix& S, Matrix& H, const std::vector<std::string>& bases, const Vector& ts, 
               const std::map<std::string, CircuitConfig>& basisDict, int numOrbitals, int matLen, bool debug = false) {
    int counter = 0;
    int counting_size = (matLen*matLen);
    //Initialize some intermediate matrices
    ComplexMatrix HMatMeasurement, HMatMeasurable, HMatExpectation, SMatMeasurement, SMatMeasurable, SMatExpectation;
    ComplexMatrix SObs, HObs;
    ComplexMatrix SMat, HMat;
    SMat.setZero();
    HMat.setZero();
    std::string nar, nal;
    for(int i = 0; i<matLen; i++){
        nal = basis_set_row[i]; //name left
        tl(-1, -1); //Default to -1 no value
        //Select approriate t, first decide if single operator or composite
        Gen_tLR(tl, nal, ts, true);
        if (nal != "R0"){
            left_er = trotter_left(nal, tl, basis_config,
                                num_spin_orbitals, num_particles, trotter_steps, atol=atol, order=order,
                                debug=False);
            left_half_h_obs = (left_er @ Ham_mat).reduce(atol=atol);
        }
        for(int j = 0; j<matLen; j++){
            counter++;
            q_ins = QuantumInstance(backend=backend, 
                                shots=int(num_shots), 
                                seed_simulator=seed,
                                seed_transpiler=seed,
                                optimization_level=opt_level,
                                noise_model=None,
                                measurement_error_mitigation_cls=mecls, cals_matrix_refresh_period=1440);

            if(print_progress)
                print("Starting {:.2f}%".format(100*(counter)/counting_size), end='\r');

            nar = basis_set_col[j] # name right
            if return_circuits:
                circ_dict[nal+'_'+nar] = {'S':[], 'H': []}
                pauli_dict[nal+'_'+nar] = {'S':0, 'H': 0}
            ## Select approriate t, first decide if single operator or composite
            if nar.count('R') == 1:
                tr = (ts[int(nar.split('R')[-1])],)
            elif nar.count('R') == 2: # composite
                nar_nopm = re.sub(r'[+-]','',nar) # remove +- signs
                nar_indlist = list(filter(None, nar_nopm.split('R'))) # remove empty strings from split
                tr = (ts[int(nar_indlist[0])], ts[int(nar_indlist[1])])
            else:
                raise ValueError('Invalid basis name for basis: {}'.format(nar))

            # NOTE: Note that for basis on the left, e.g., <+R1-R2| in <+R1-R2|H|+R5-R5>, t's are reversed
            # here, so tl=(t2,t1) in the above example, BUT NOT the BASIS NAME, the name is still '+R1-R2'
            # the reverse of the name is in base_circuit_code_separate(), 
            # and adding - sign to tl is in base_circuit_code_single()
            if nal == nar and nal == 'R0': # Diagonal Entries
                ## Compute S matrix Entry
                S_mat[i,j] = 1
                ## Compute H matrix Entry
                sampler_h = CircuitSampler(q_ins)
                h_mat_measurement = StateFn(Ham_mat).adjoint()
                h_mat_measurable_expression = h_mat_measurement  @ CircuitStateFn(init_state)
                h_mat_trotterized_op = PauliTrotterEvolution(trotter_mode=Suzuki(reps=1, order=order)).convert(h_mat_measurable_expression)
                h_mat_expectation = PauliExpectation(group_paulis=True).convert(h_mat_trotterized_op)
                h_mat_sampled_exp_op = sampler_h.convert(h_mat_expectation)
                H_mat[i,j] = h_mat_sampled_exp_op.eval()
                if print_entry:
                    print(nal,',', nar,':', S_mat[i,j],',', H_mat[i,j])
                # Save the circuit
                if return_circuits:
                    pauli_dict[nal+'_'+nar]['H'] = len(Ham_mat)
                    h_circ_ops = list(sampler_h._circuit_ops_cache.values())
                    for hc in h_circ_ops:
                        circ_dict[nal+'_'+nar]['H'].append(hc.to_circuit())
            else:
                if nar == 'R0':
                    s_obs = left_er
                    h_obs = left_half_h_obs
                elif nal == 'R0':
                    right_er = trotter_right(nar, tr, basis_config,
                                    num_spin_orbitals, num_particles, trotter_steps, atol=atol, order=order,
                                    debug=False)  
                    s_obs = right_er
                    h_obs = (Ham_mat @ right_er).reduce(atol=atol)
                else:
                    right_er = trotter_right(nar, tr, basis_config,
                                    num_spin_orbitals, num_particles, trotter_steps, atol=atol, order=order,
                                    debug=False)  
                    s_obs = (left_er @ right_er).reduce(atol=atol)
                    h_obs = (left_half_h_obs @ right_er).reduce(atol=atol)
                ## Compute S matrix Entry
                sampler_s = CircuitSampler(q_ins)
                s_mat_measurement = StateFn(s_obs).adjoint()
                s_mat_measurable_expression = s_mat_measurement  @ CircuitStateFn(init_state)
                s_mat_expectation = PauliExpectation(group_paulis=True).convert(s_mat_measurable_expression)
                s_mat_sampled_exp_op = sampler_s.convert(s_mat_expectation)
                S_mat[i,j] = s_mat_sampled_exp_op.eval()
                if symmetric:
                    S_mat[j,i] = S_mat[i,j].conjugate()
                ## Compute H matrix Entry
                sampler_h = CircuitSampler(q_ins)
                h_mat_measurement = StateFn(h_obs).adjoint()
                h_mat_measurable_expression = h_mat_measurement  @ CircuitStateFn(init_state)
                h_mat_expectation = PauliExpectation(group_paulis=True).convert(h_mat_measurable_expression)
                h_mat_sampled_exp_op = sampler_h.convert(h_mat_expectation)
                H_mat[i,j] = h_mat_sampled_exp_op.eval() # double
                if symmetric:
                    H_mat[j,i] = H_mat[i,j].conjugate()
                if print_entry:
                    print(nal,',', nar,':', S_mat[i,j],',', H_mat[i,j])
                # ## Save the circuit
                if return_circuits:
                    pauli_dict[nal+'_'+nar]['S'] = len(s_obs)
                    pauli_dict[nal+'_'+nar]['H'] = len(h_obs)
                    s_circ_ops = list(sampler_s._circuit_ops_cache.values())
                    h_circ_ops = list(sampler_h._circuit_ops_cache.values())
                    for sc in s_circ_ops:
                        circ_dict[nal+'_'+nar]['S'].append(sc.to_circuit())
                    for hc in h_circ_ops:
                        circ_dict[nal+'_'+nar]['H'].append(hc.to_circuit())
    
    }
    S = SMat;
    H = HMat;
    //if return_circuits:
        //return S_npmat, H_npmat, circ_dict, pauli_dict
}