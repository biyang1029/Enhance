#include "DataConfig.h"
#include "SimulationController.h"
#include <algorithm>
#include <cctype>
#include <iostream>


int main() {
    // Optional: choose refrigerant & key HP params at startup (press Enter to keep default)
    try {
        std::string in;
        std::cout << "Use REFPROP for refrigerant properties? [Y/n] (default Y): ";
        std::getline(std::cin, in);
        if (!in.empty()) {
            char c = static_cast<char>(std::tolower(static_cast<unsigned char>(in.front())));
            CFG.hp.use_refprop = (c != 'n');
        }
        if (CFG.hp.use_refprop) {
            std::cout << "REFPROP directory (default " << CFG.hp.refprop_dll_dir << "): ";
            std::getline(std::cin, in);
            if (!in.empty()) CFG.hp.refprop_dll_dir = in;
        } else {
            CFG.hp.refprop_dll_dir.clear();
        }

        std::cout << "Select refrigerant [R134a/R410A/R32/R290] (default R134a): ";
        std::getline(std::cin, in);
        if (!in.empty()) CFG.hp.fluid = in;
        std::cout << "Isentropic efficiency (default 0.70): ";
        std::getline(std::cin, in);
        if (!in.empty()) CFG.hp.eta_isentropic = std::max(0.3, std::min(0.9, std::stod(in)));
        std::cout << "Superheat K (default 3): ";
        std::getline(std::cin, in);
        if (!in.empty()) CFG.hp.superheat_K = std::max(0.0, std::stod(in));
        std::cout << "Subcool K (default 5): ";
        std::getline(std::cin, in);
        if (!in.empty()) CFG.hp.subcool_K = std::max(0.0, std::stod(in));
    } catch (...) { /* keep defaults */ }

    SimulationController sim(CFG);
	sim.run("results.csv");
	return 0;
}
