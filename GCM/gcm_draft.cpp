#include "gcm_header.hpp"

int main() {

    std::vector<std::string> bases = {
        "R0",
        "+R1", "-R1",
        "+R2", "-R2",
        "+R3", "-R3",
        "+R4", "-R4",
        "+R3+R5", "+R4+R6",
        "+R2+R1", "+R2-R1", "-R2+R1", "-R2-R1"
    };
    std::map<std::string, CircuitConfig> basisDict = {
        {"R1", {{"5^ 2", 1}, {"7^ 4", -1}, {"2^ 5", -1}, {"4^ 7", 1}}},
        {"R2", {{"3^ 0", 1}, {"6^ 1", -1}, {"0^ 3", -1}, {"1^ 6", 1}}},
        {"R3", {{"1^ 0", 1}, {"2^ 3", -1}, {"0^ 1", -1}, {"3^ 2", 1}}},
        {"R4", {{"6^ 3", 1}, {"7^ 5", -1}, {"3^ 6", -1}, {"5^ 7", 1}}},
        {"R5", {{"0^ 4", 1}, {"2^ 6", -1}, {"4^ 0", -1}, {"6^ 2", 1}}},
        {"R6", {{"1^ 3", 1}, {"3^ 5", -1}, {"5^ 1", -1}, {"7^ 4", 1}}}
    };
    

    std::cout << "Test";

    //Generate t's
    uint64_t seed = 1234; //Seed for the PCG64 generator
    Vector ts(7); //Generate random numbers w/ PCG64 and scale them to be between 0 and 100
    ts[0] = 1.0; //Inserting 1 for R0
    for (int i = 1; i < ts.size(); i++) {
        ts[i] = generateRandomNumber(seed);
    }
    

    Matrix Hamiltonian; // Hamiltonian matrix
    Matrix HFState;     // Hartree-Fock state matrix, used as the initial input
  
    int numOrbitals = 8; // Define number of orbitals in the system
    
    //Step 1: Perform Jordan-Wigner transformation
    //Step 2: Generate unitaries
    //Step 3: Trotterize unitaries into Pauli strings
    //Step 4-8: Compute matrix elements classically
    int matLen = bases.size();
    int matSize = matLen * matLen;
    //Define matrices of complex numbers
    Matrix Sclass(matLen, matLen);
    Matrix Hclass(matLen, matLen);
    //Initialize matrices with zeros
    Sclass.setZero();
    Hclass.setZero();
    //Compute function fills in the S and H matrices
    ComputeSH(HFState, Hamiltonian, Sclass, Hclass, bases, ts, basisDict, numOrbitals, matLen, true);
    std::cout << Sclass;
    std::cout << Hclass;

    //Steps 9 & 10: Evaluate S and H quantumly
    Matrix Squant(matLen, matLen);
    Matrix Hquant(matLen, matLen);
    //Initialize matrices with zeros
    Squant.setZero();
    Hquant.setZero();
    QuantumSH(HFState, Hamiltonian, Sclass, Hclass, bases, ts, basisDict, numOrbitals, matLen, true);

    //Step 13: Solve the general eigenvalue problem classically
    //Matrix overlapMatrix;
    //auto [eigenvalues, eigenvectors] = SolveEigenvalueProblem(evaluatedMatrixElements, overlapMatrix);

    //Output the results
    //std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;
    //std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
    return 0;
}