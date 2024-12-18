#include <iostream>

#include "xacc.hpp"
#include "nwq_accelerator.hpp"

int main(int argc, char **argv)
{
  xacc::Initialize(argc, argv);
  // Accelerator:

  auto nwq_acc = std::make_shared<xacc::quantum::NWQAccelerator>();

  if (argc == 2)
    nwq_acc->updateConfiguration({std::make_pair("backend", std::string(argv[1]))});

  auto qpp_acc = xacc::getAccelerator("qpp"); //, {std::make_pair("shots", 4096)}

  xacc::qasm(R"(
  .compiler xasm
  .circuit bell
  .qbit q
  H(q[0]);
  CX(q[0], q[1]);
  Measure(q[0]);
  Measure(q[1]);
  )");
  auto bell = xacc::getCompiled("bell");

  // Allocate some qubits and execute
  auto buffer_nwq = xacc::qalloc(2);
  nwq_acc->execute(buffer_nwq, bell);
  std::cout << "NWQ (Counts): "
            << "\n";
  buffer_nwq->print();

  // Allocate some qubits and execute
  buffer_nwq->resetBuffer();
  nwq_acc->updateConfiguration({std::make_pair("vqe_mode", true)});
  nwq_acc->execute(buffer_nwq, bell);
  std::cout << "NWQ (EXP-Z): "
            << "\n";
  buffer_nwq->print();

  auto buffer_qpp = xacc::qalloc(2);
  qpp_acc->execute(buffer_qpp, bell);
  std::cout << "QPP: "
            << "\n";
  buffer_qpp->print();

  xacc::Finalize();
  return 0;
}