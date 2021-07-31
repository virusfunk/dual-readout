#include "fastjetInterface.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include <algorithm>
#include <cmath>

fastjetInterface::fastjetInterface() {
  fJets = new std::vector<fastjetData>(0);
  fJetBase = new std::vector<fastjetDataBase>(0);
}

fastjetInterface::~fastjetInterface() {
  if (fJets) delete fJets;
  if (fJetBase) delete fJetBase;
}

fastjetInterface::fastjetDataBase::fastjetDataBase(fastjet::PseudoJet& jet) {
  E = jet.E();
  px = jet.px();
  py = jet.py();
  pz = jet.pz();
}

fastjetInterface::fastjetData::fastjetData(fastjet::PseudoJet& jet) {
  E = jet.E();
  px = jet.px();
  py = jet.py();
  pz = jet.pz();
  phi = jet.phi();
  phi_std = jet.phi_std();
  rap = jet.rap();
  eta = jet.eta();
  pt = jet.pt();
  m = jet.m();
  mt = jet.mt();
  hasAssociatedCS = jet.has_associated_cs();
  validCS = jet.has_valid_cs();
  hasConstituents = jet.has_constituents();
  nConstituents = jet.constituents().size();

  fastjet::PseudoJet childPJ;
  hasChild = jet.has_child(childPJ);
  child = fastjetDataBase(childPJ);
}

void fastjetInterface::init(TTree* treeIn, std::string branchname) {
  treeIn->Branch(branchname.c_str(), &fJets);
}

void fastjetInterface::writeJets(std::vector<fastjet::PseudoJet> jets) {
  fJets->clear();
  fJets->reserve(jets.size());

  for (auto jet = jets.begin(); jet != jets.end(); ++jet) {
    fJets->push_back(fastjetData(*jet));
  }
}

void fastjetInterface::MywriteJets(std::vector<fastjet::PseudoJet> jets, std::vector<fastjet::PseudoJet> constituents[], std::vector<long long int> constituents_fibernum[]) {
    fJets->clear();
    fJets->reserve(jets.size());

    for (auto jet = jets.begin(); jet != jets.end(); ++jet) {
        fJets->push_back(fastjetData(*jet));
    }

    fjets_constituents = new std::vector<fastjetData>[num_jets];
    for (int i = 0; i < num_jets; i++) {
        for (int j = 0; j < constituents[i].size(); j++) {
            fjets_constituents[i].push_back(fastjetData(constituents[i].at(j)));
        }
    }

    fjets_constituents_fibernum = new std::vector<long long int>[num_jets];
    for (int i = 0; i < num_jets; i++) {
        for (int j = 0; j < constituents_fibernum[i].size(); j++) {
            fjets_constituents_fibernum[i].push_back(constituents_fibernum[i].at(j));
        }
    }
}


void fastjetInterface::runFastjet(const std::vector<fastjet::PseudoJet>& input) {
  // FastJet
  double dR = 0.8;
  fastjet::JetDefinition jetDef(fastjet::ee_genkt_algorithm,dR,-1);

  // Run Fastjet algorithm
  std::vector<fastjet::PseudoJet> inclusiveJets, sortedJets;
  fastjet::ClusterSequence clustSeq(input, jetDef);

  // Extract inclusive jets sorted by pT
  inclusiveJets = clustSeq.inclusive_jets();
  sortedJets    = fastjet::sorted_by_pt(inclusiveJets);

  writeJets(sortedJets);
}

void fastjetInterface::MyrunFastjet(const std::vector<fastjet::PseudoJet>& input, const std::vector<long long int>& fibernum) {
    // FastJet
    double dR = 0.8;
    fastjet::JetDefinition jetDef(fastjet::ee_genkt_algorithm, dR, -1);

    // copy input
    std::vector<fastjet::PseudoJet> editable_input = input;

    // put index
    for (int i = 0; i < input.size(); i++) {
        editable_input.at(i).set_user_index(i);
    }

    // Run Fastjet algorithm
    std::vector<fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(editable_input, jetDef);

    // Extract inclusive jets sorted by pT
    inclusiveJets = clustSeq.inclusive_jets();
    sortedJets = fastjet::sorted_by_pt(inclusiveJets);

    num_jets = sortedJets.size();
    std::vector<fastjet::PseudoJet> constituents[num_jets];
    for (int i = 0; i < num_jets; i++) {
        constituents[i] = clustSeq.constituents( sortedJets.at(i) );
    }

    std::vector<long long int> constituents_fibernum[num_jets];
    for (int i = 0; i < num_jets; i++) {
        for (auto constituent : constituents[i]) {
            int temp_index = constituent.user_index();
            constituents_fibernum[i].push_back(fibernum.at(temp_index));
        }
    }
    // exception
    if (editable_input.size() != fibernum.size()) {
        printf("the size of pseudojet is different from fibernum's size\n");
        printf("size of input: %d\n", input.size());
        printf("size of editable input: %d\n", editable_input.size());
        printf("size of fibernum input: %d\n", fibernum.size());
        std::abort();
    }
    for (int i = 0; i < num_jets; i++) {
        if (constituents[i].size() != constituents_fibernum[i].size()) {
            printf("the size of constituents is different from fibernum's size\n");
            std::abort();
        }
    }
    for (int i = 0; i < num_jets; i++) {
        for (auto constituent : constituents[i]) {
            if (!sortedJets.at(i).contains(constituent)) {
                printf("cluster does not contain a constituent\n");
                std::abort();
            }
        }
    }

    MywriteJets(sortedJets, constituents, constituents_fibernum);
}


void fastjetInterface::set(TTree* treeIn, std::string branchname) {
  treeIn->SetBranchAddress(branchname.c_str(),&fJets);
}

void fastjetInterface::read(std::vector<fastjetData>& jets) {
  jets = *fJets;
}

void fastjetInterface::Myread(std::vector<fastjetData>& jets, std::vector<fastjetData>** constituents, int* num, std::vector<long long int>** constituents_for_fibernum) {
    jets = *fJets;
    *constituents = fjets_constituents;
    *num = num_jets;
    *constituents_for_fibernum = fjets_constituents_fibernum;
}
