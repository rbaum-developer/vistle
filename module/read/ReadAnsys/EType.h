/* This file is part of COVISE.

   You can use it under the terms of the GNU Lesser General Public License
   version 2.1 or later, see lgpl-2.1.txt.

 * License: LGPL 2+ */

#ifndef _ETYPE_H_
#define _ETYPE_H_

struct EType {
    int id_; // ITEM 1: Referenz Nummer des Elements
    int routine_; // ITEM 2:Element routine Nummer
    int keyops_[12]; // ITEM 3-14: nicht näher definierte Optionen
    int dofpernode_; // ITEM 34: Anzahl der DOFs pro node
    int nodes_; // ITEM 61: Anzahl der Knoten für dieses Element
    int nodeforce_; // ITEM 63: Anzahl der Knoten mit Kräften
    int nodestress_; // ITEM 94: Anzahl der Knoten mit Schubspannungen
};
#endif
