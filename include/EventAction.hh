//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// Phys. Med. Biol. 63(10) (2018) 105014-12pp
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file PhysicsList.hh
/// \brief Definition of the PhysicsList class

#ifndef EventAction_hh
#define EventAction_hh 1

#include "G4DNAChemistryManager.hh"
#include "G4UserEventAction.hh"
#include "G4Version.hh"

class EventAction : public G4UserEventAction
{
  public:
    void BeginOfEventAction(const G4Event* event) override
    {
#if G4VERSION_NUMBER >= 1140
      if (G4DNAChemistryManager::GetInstanceIfExists() != nullptr)
        G4DNAChemistryManager::Instance()->BeginOfEventAction(event);
#endif
    }
    void EndOfEventAction(const G4Event* event) override
    {
#if G4VERSION_NUMBER >= 1140
      if (G4DNAChemistryManager::GetInstanceIfExists() != nullptr)
        G4DNAChemistryManager::Instance()->EndOfEventAction(event);
#endif
    }
};

#endif
