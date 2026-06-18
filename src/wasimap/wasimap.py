"""
This file contains the main WaSiMap logic

Pedro Manuel Guillem Gloria (2026)
Pisabarro Group
Structural Bioinformatics, BIOTEC
TU Dresden
"""

from __future__ import annotations
import mdtraj as md
import numpy as np
from multiprocessing import Process
from multiprocessing import Manager 
from pathlib import Path
from datetime import datetime

from importlib.resources import files
from pathlib import Path
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from functools import partial
import threading, webbrowser, json, shutil, time, sys, requests, os

#MDTraj compatible formats (to March 2026)
TRAJECTORY_EXTS = {".xtc", ".trr", ".dcd", ".nc", ".h5", ".hdf5", ".lh5", ".binpos"}
TOPOLOGY_EXTS   = {".pdb", ".gro", ".prmtop", ".parm7", ".psf", ".top"}
H5_LIKE         = {".h5", ".hdf5", ".lh5"}

HBOND_ATOMS = {"O", "N", "S"} #Atoms usually mediating hbonding.. Sulfur is weak.. all others are negligible

#For future use - Handle solvents other than just water
SOLVENTS = {
    # --- Water ---
    "HOH", "WAT", "SOL", "TIP3", "TIP3P", "TIP4P", "TIP4", "SPC", "SPCE", "OPC3", "OPC", "TIP5"

    # --- Simple alcohols ---
    "MET", "MEOH", "METHANOL",
    "ETH", "EOH", "ETOH", "ETHANOL",
    "PRO", "IPR", "IPRO", "ISOPROPANOL",
    "BUT", "BUOH", "TBU", "TBA", "TERTBUTANOL",

    # --- Halogenated solvents ---
    "CL3", "CHCL3", "CHL", "CHCL", "CHLOROFORM",
    "DCM", "CH2CL2", "DICHLOROMETHANE",
    "DCE", "CH2CLCH2CL", "DICHLOROETHANE",
    "TCE", "TRICHLOROETHANE",
    "CTC", "CCl4", "CARBON_TETRACHLORIDE",

    # --- Aprotic polar solvents ---
    "DMS", "DMSO",
    "DMF",
    "DMA", "DMAC",
    "ACN", "CH3CN", "ACETONITRILE",
    "ACE", "ACETONE",
    "NMP",

    # --- Ethers ---
    "THF",
    "DIOX", "DIOXANE",
    "DEE", "DIETHYLETHER", "ETHER",
    "MTBE",

    # --- Aromatic solvents ---
    "BEN", "BENZENE",
    "TOL", "TOLUENE",
    "XYL", "XYLENE",

    # --- Alkanes ---
    "HEX", "HEXANE",
    "HEP", "HEPTANE",
    "OCT", "OCTANE",
    "NON", "NONANE",
    "DEC", "DECANE",
    "CYC", "CYCLOHEXANE",

    # --- Misc common ---
    "GLY", "GLYCOL", "ETHYLENEGLYCOL",
    "PG", "PROPYLENEGLYCOL",
    "UREA",
    "FORM", "FORMAMIDE",

    # --- Cryo / crystallography additives ---
    "PEG", "PEG400", "PEG3350",
    "MPD",
    "GOL", "GLC", "GLYCEROL",
    "SO4", "PO4",  # sometimes behave like solvent-like species
}


# This class extends SimpleHTTPRequestHandler, it 
# overrides end_headers to force no cache, and log_message() 
# to dump web server events to stdout
class NoCacheHandler(SimpleHTTPRequestHandler):

    def end_headers(self):
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def log_message(self, format, *args):
        print("[wasimap]", format % args, file=sys.stdout, flush=True)


#/////////////////////////////////////////////////////////////////////////////
#
# Calculates euclidean distances from every water to a given heavy-atom, in every frame.
# I do this in a function so we can enqueue with a process pool in parallel
#
# id:         atomID of heavy-atom
# anchordict: the result array
# waters:     IDs of every water oxygen in the trajectory
# trajxyz:    XYZ numpy array of traj object (passing the entire traj would eat up the memory)
# around_nanometers: cuttoff distance
#
def processAtom(id, anchordict, waters, trajxyz, around_nanometers):
    #for id in involved_ids:
    print(f"Processing AtomID {id}")
    wadict = {}
    for water in waters:
        #This is the magic.. we use numpy to compute euclideans from all waters to a heavy-atom, for all frames, in a flash :)
        framebyframe = np.sqrt(np.sum((trajxyz[:, water, :] - trajxyz[:, id, :])**2, axis=1))
        framenum = 0
        selectedframes = [] #Stores coordinates of frames
        for euclidean in framebyframe: #Iterate euclidean distance of waters to IDS in the matrix
            if euclidean <= around_nanometers:
                    if wadict.get(water) == None:
                        wadict[water] = []
                    selectedframes.append(framenum)
            framenum=framenum+1
        wadict[water] = selectedframes
        if wadict.get(water) != None:
            if len(wadict[water]) < 20: #Ditch meaningless water.. we assume if the water occupies a site less than 50 frames, is not representative
                wadict.pop(water)
    anchordict[id] = wadict
    #print(anchordict)
    print(f"##############FINISHED HEAVY-ATOM {id}###############")


class WaterMapper:

    around_nanometers = 0
    relevance         = 0 
    n_frames          = 0

    
    #Constructor... Watch out.. distance threshold comes in nanometers, not angstroms
    def __init__(self, distance_threshold=0.35, persistence=5, gui=False, onlygui=False, inputs=None, testdata=False, output_folder="./wasimap_outputs",):
        #Assign Variables
        self.around_nanometers  = distance_threshold #In nanometers
        self.relevance          = persistence/100 #in percentage (int)
        self.path               = "./"
        self.gui                = gui
        self.onlygui            = onlygui
        self.testdata            = testdata
        self.output_folder      = output_folder
        self.inputs             = inputs or "auto"
        self.md_trajectories    = {}

    def run(self) -> None:
        print("Running WaSiMap")
        print(f"  distance_cutoff={self.around_nanometers} nm")
        print(f"  persistence_threshold={self.relevance*100}%")
        print(f"  gui={self.gui}")
        print(f"  onlygui={self.onlygui}")
        print(f"  inputs={self.inputs}")
        print(f"  testdata={self.testdata}")
        print(f"")
        print(f"Structural Bioinformatics Group. BIOTEC. TU Dresden 2026")
        print(f"If you find our work useful, please cite us")

        #If testdata is to be downloaded
        #Test Data from public Zenodo repository DOI https://doi.org/10.5281/zenodo.18984212
        if(self.testdata):
          d = Path(self.output_folder)
          print("")
          print("Test data resides at Zenodo, with DOI https://doi.org/10.5281/zenodo.18984212")
          print("WaSiMap will download 320Mb of test MD trajectory data to the current folder")
          print("")
          answer = input("Would you like to continue? [y/N]: ").strip().lower()
          if answer not in ("y", "yes"):
                print("Aborted by user.")
                raise SystemExit(0)
          else:
            #URL for the 4 files (2 H5 trajectories and 1 NC trajectory, with prmtop topology)
            urls = [
                "https://zenodo.org/records/18984212/files/sim1.h5?download=1",
                "https://zenodo.org/records/18984212/files/sim2.h5?download=1",
                "https://zenodo.org/records/18984212/files/sim3.nc?download=1",
                "https://zenodo.org/records/18984212/files/sim3.prmtop?download=1",
            ]

            for url in urls:
                filename = os.path.basename(url.split("?")[0])
                print(f"Downloading {filename}... from {url}")

                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    with open(filename, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)

            print("All downloads completed") 
            print("")
            print("Now you can run 'wasimap --gui'")      
            print("")
            exit(0)




        #If this is a GUI only run.. execute server
        if(self.onlygui):
          d = Path(self.output_folder)
          if d.exists():
            a = Path(f"{self.output_folder}/wasimap.json")
            if a.exists():
                self.serveGUI("0.0.0.0",8080)
            else:
                print("Output folder doesn't contain results")
          else:
            print("Output folder doesn't exist")
          print("Program ended")
          exit(0)


        #Start looking for MD trajectories on the local directory
        if self.inputs == "auto":
          print(f"-----------------------------------------------")
          print(f"Scanning local directory for MD trajectory data")
          print(f"-----------------------------------------------")
          print()
          self.md_trajectories = self.scanForTrajs(".") #Scan auto
        else:
          self.md_trajectories = self.resolve_md_pair(self.inputs) #check traj and top exist with given name



        #Enumerate & confirm if user wants to proceed
        print("Found the following trajectory / topology pairs:\n")
        for key, value in self.md_trajectories.items():
            if value:
                print(f"{key}  -->  {value}")
            else:
                print(f"{key}  -->  [embedded topology or not required]")
        print()
        answer = input("Proceed with these files? [y/N]: ").strip().lower()
        if answer not in ("y", "yes"):
            print("Aborted by user.")
            raise SystemExit(0)
        


        #Ask user for desired interfacial positions to watch (common to all trajectories)
        positions_raw = input("Manual residue number or range to include (example: 1-6 or 1,2,3,N) (Press Enter to auto-detect): ").strip()
        positions = [
            i
            for part in positions_raw.split(",")
                for i in (
                    range(int(part.split("-", 1)[0]), int(part.split("-", 1)[1]) + 1)
                    if "-" in part
                    else [int(part)]
                )
        ] if positions_raw else []


        #Ask user if the entire complex is in a single chain
        last_residue_raw = input(
            "If known, provide the last residue number of the first molecule (Press Enter to auto-detect): "
        ).strip()

        if last_residue_raw == "":
            last_residue = None
        else:
            try:
                last_residue = int(last_residue_raw)
            except ValueError:
                raise SystemExit("Residue number must be an integer.")


        #Prepare exit dictionary
        self.prepare_outdir()

        #Create output dictionary structure
        bundle = {
            "simulations": {},   # dict of many sim results
            "config": {"distance_cutoff": self.around_nanometers, "persistence_threshold": self.relevance*100},        # optional: cutoff, thresholds, etc.
            "meta": {"generator": "WaSiMap","timestamp": datetime.utcnow().isoformat()},          # optional: timestamps, versions, etc.
        }

        #FOR EACH TRAJECTORY FOUND, START MAIN ROUTINE
        for traj, topo in self.md_trajectories.items():
           print(f"-----------BEGIN MD SIMULATION ANALYSIS FOR {traj}----------------")
           bundle["simulations"][traj] = self.findWetSpots(traj, topo, positions, last_residue) 
        
        #Save JSON results on output folder
        with open(f"{self.output_folder}/wasimap.json", "w+") as f:
            json.dump(bundle, f, indent=2)
        
        #Generate HTML5 website with results
        #Copy html viewer from resource into users toplevel directory
        if(self.gui):
            html = files("wasimap").joinpath("wasimap.html").read_text(encoding="utf-8")
            Path(f"{self.output_folder}/wasimap.html").write_text(html, encoding="utf-8")
            self.serveGUI("0.0.0.0",8080)
        
        print()
        print()
        print("If you used WaSiMap in your research, please cite our work :-)")
        
        exit(0)










    #Start a small web server to open GUI properly. Page must source local files
    def serveGUI(self, host="0.0.0.0", port=8080):
        server = ThreadingHTTPServer((host, port), partial(NoCacheHandler, directory=self.output_folder))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        webbrowser.open(f"http://localhost:{port}/wasimap.html")
        
        print(f"Server running at http://{host}:{port}")
        print("Press CTRL+C to stop")

        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("\nStopping server...")
            server.shutdown()
            server.server_close()
        return server



    #Create output directory for files
    def prepare_outdir(self):
        d = Path.cwd() / self.output_folder
        if d.exists():
            shutil.rmtree(d)
        d.mkdir()
        return d



    #Scan for trajectory/topology combos in the current directory
    def scanForTrajs(self, d="."):
        d = Path(d)
        files = [f for f in d.iterdir() if f.is_file()]
        topo_by_stem = {f.stem: f.name for f in files if f.suffix.lower() in TOPOLOGY_EXTS}
        trajs = [f for f in files if f.suffix.lower() in TRAJECTORY_EXTS]

        if not trajs:
            raise FileNotFoundError(f"No trajectory files found in current directory")

        pairs = {}
        missing = []

        for t in trajs:
            ext = t.suffix.lower()
            if ext in H5_LIKE:
                pairs[t.name] = ""
            else:
                top = topo_by_stem.get(t.stem)
                if top:
                    pairs[t.name] = top
                else:
                    missing.append(t.name)

        if not pairs and not missing:
            raise FileNotFoundError(f"No usable trajectory/topology pairs found in: {d}")

        if missing:
            raise SystemExit("Topologies are needed for trajectories: " + ", ".join(sorted(missing)))

        return pairs



    #Check if given input name exists as trajectory and topology. Return a dict with both names if it does
    def resolve_md_pair(stem: str, directory: str | Path = "."):
        d = Path(directory)
        pair = {}
        traj = next((d / f"{stem}{ext}" for ext in TRAJECTORY_EXTS if (d / f"{stem}{ext}").exists()), None)
        top  = next((d / f"{stem}{ext}" for ext in TOPOLOGY_EXTS   if (d / f"{stem}{ext}").exists()), None)
        if not traj:
            raise FileNotFoundError(f"No trajectory found for '{stem}'")
        if not top:
            raise FileNotFoundError(f"No topology found for '{stem}'")
        
        pair[traj] = top
        return pair
    



   #-------------------------------------------------------------------------------------------------------
   #Main Mapping method starts here
    def findWetSpots(self, trajectory, topology, userpositions = [], last_residue = ""):       
            print(f"Processing Trajectory {trajectory} and topology {topology}")
            print(f"----------------------------------")

            #load it - If H5, single file - if other types, load topology too (thank you MDTRAJ!!!)
            if Path(trajectory).suffix.lower() in H5_LIKE:
             traj = md.load(f"./{trajectory}")
            else:
             traj = md.load(f"./{trajectory}", top=f"./{topology}")
            
            print(f"Imaging trajectory to single periodic box")        
            # 1. Image molecules (fix broken molecules due to PBC)
            traj.image_molecules(inplace=True)
            # 2. Center coordinates (move system to origin, solvent included)
            traj.center_coordinates()
            # 3. Align entire system to the first frame (solvent + protein)
            traj.superpose(traj, frame=0)

            around_nanometers = self.around_nanometers #cuttof distance from water to heavy atom
            nframes           = traj.n_frames #number of frames in the trajectory
            relevant          = nframes*self.relevance #minimum residence frames to be considered a relevant water
            involved_ids      = []
            
            #Extract first frame of traj, and save 
            traj.atom_slice(traj.topology.select("not resname HOH and not resname WAT and not resname SOL"))[0].save_pdb(f"{self.output_folder}/wasimap_{trajectory}.pdb")
            #GEt IDS of atoms involved in the interface
            
            #Add atoms of user defined positions
            print(f"Adding user defined positions (if any)")
            for pos in userpositions:
              res = int(pos)-1
              #Select heavy atoms of position
              ap = traj.topology.select(f"(resid {res}) and ((symbol == 'O') or (symbol == 'N') or (symbol == 'S'))")
              for addi in ap:
                print(f"Adding user-defined position {addi}")
                involved_ids.append(addi)

            if len(involved_ids) > 0:
               print(f"Added {len(involved_ids)} heavy atoms to the list of atom IDs to process")
               #Include intermolecular atoms that are making hbonds in the interface

            #Try to guess the last residue of first chain
            #The logic here is that two chains may exist in the simulation
            #This is only necesary if the user didn't provide the position of the last residue of a simulation with all residues incrementally in a single chain
            if last_residue is None:
                last_residue = self.find_first_molecule_end_by_peptide_break(traj, 2.0)
                print(f"First molecule ends at residue.index {last_residue} or position {int(last_residue)+1}")

            # BAKER HUBBARD HBOND FORMING ATOMS.. WE ARE INTERESTED IN THESE
            print(f"Finding atoms that form intermolecular Hbonds")
            hbonds = md.baker_hubbard(traj, periodic=True, sidechain_only=False, freq=self.relevance)

            print(f"Found {len(hbonds)} hbonds with Baker-Hubbard method.. scanning heavy-atoms")
            for hbond in hbonds:
                #BUILD ATOM COLLECTION
                atom1 = traj.topology.atom(hbond[0])
                atom2 = traj.topology.atom(hbond[2])
                if atom1.residue.index+1 < last_residue and atom2.residue.index+1 < last_residue: 
                    continue
                elif atom1.residue.index+1 > last_residue and atom2.residue.index+1 > last_residue:
                    continue
                else:
                    involved_ids.append(hbond[0])
                    involved_ids.append(hbond[2])
                    print(f"Adding atoms {hbond[0]} and {hbond[2]}")

            print(f"Added heavy atoms that form inter-molecular Hbonds")

            if len(involved_ids) == 0:
                print(f"No heavy atoms detected for analysis :/.. Are your molecules too far appart?. Provide a residue list manually")
                res = {'filename': trajectory, 'anchor_contacts' : {}, 'important_waters': {}}
                return res 


            ##processAtom was here.. I moved outside of the class so windows can spawn it (instead of forking)
            ##Linux seems to be fine with forking


            #Launch Processes in Parallel
            manager = Manager() #Create a manager processing pool
            involved_ids = list(set(involved_ids)) #Remove duplicates from atom list
            print(f"Involved AtomIDs {involved_ids}")

            
            #We could build selectors here to work with any type of solvent (chloroform, etc)
            #And get atomic IDs with hbonding potential, including water oxygens for dynamics with 2 or more solvents + water

            #mask = (
            #    "(" + " or ".join(f"resname {r}" for r in sorted(SOLVENTS)) + ") and "
            #    "(" + " or ".join(f"symbol {e}" for e in sorted(HBOND_ATOMS)) + ")"
            #)

            #GET ids of Oxygens in water atoms
            mask   = "(water) and (symbol == 'O')" #only water for now
            waters = traj.topology.select(mask) #only water for now
            
            waterdict = {}
            anchordict = manager.dict() #Shared proxy object.. we need child forks to propagate back
            jobs    = [] #Process pool
            
            #Create a process per AtomID
            for id in involved_ids:
                anchordict[id] = {}
                p = Process(target=processAtom, args=(id, anchordict, waters, traj.xyz, around_nanometers))
                jobs.append(p)
                p.start()

            #Launch parallel jobs    
            for proc in jobs:
                proc.join()
                        
            #print("this is what came out")
            #print(anchordict)
            #Define a list of important water ids
            importantwaters = {} #Stores all important waterIDs, frame arrays near, and residence percentage
            anchorcontacts  = {} #Stores anchorpoint and for each the waterids that crossed it
            watersites = {} #Stores the average coordinates of the watersites, and the waters involved
            watersiteid = 1
            for anchor, waters in anchordict.items():
                anchor = str(anchor)
                owat = {k: v for k, v in sorted(waters.items(), key=lambda item: item[1], reverse=True)}
                print(f"{len(waters)} waters crossed near anchor-atom {anchor}")
                aguas = []
                for atom, frames in owat.items(): #atom is the water atom ID, frames is an array of frames where water was near anchor
                    if len(frames) >= relevant:
                        print(f"Appending {atom} to list of importants ({len(frames)} frames)")
                        
                        if importantwaters.get(str(atom)) == None:
                            importantwaters[str(atom)] = {}
                        
                        importantwaters[str(atom)]['residue'] = str(traj.topology.atom(atom)) #The water residue (not atom id)
                        importantwaters[str(atom)]['frames']  = frames #array of frame indexes where water was near
                        #PEDRO 2025 - ADD OCCUPANCY/RESIDENCE 
                        atom_coords = traj.xyz[frames, atom, :]  # Shape: (selected_frames, water_atomid) # Gets an np array of 3d coords of water at selected frames
                        # Compute the average position (x, y, z)
                        average_position = np.mean(atom_coords, axis=0)#get the mean 3D position of this water (defines a water site)
                        average_position *= 10 #convert to Angstroms (units are in nanometers)
                        importantwaters[str(atom)]['residence_percentage']  = round((len(frames)/nframes)*100) #Get the residence time of this water
                        importantwaters[str(atom)]['wetspot']  = average_position.tolist() #Numpy array
                        #x_avg, y_avg, z_avg = average_position 
                        #print(f"X  : {x_avg}, Y: {y_avg}, Z: {z_avg}")
                        if anchorcontacts.get(anchor) == None:
                            anchorcontacts[anchor] = {}
                        
                        anchorcontacts[anchor]['residue'] = str(traj.topology.atom(int(anchor))) #anchorpoint residue in RESPOS-ATOM format
                        aguas.append(int(atom)) 
                if len(aguas) > 0:
                    
                    #BUILD WATERSITES
                    coords = []
                    for agua in aguas:
                           coords.append(importantwaters[str(agua)]['wetspot'])
                    centroid = np.mean(np.array(coords), axis=0)
                    centroid = centroid.tolist()
                    #BUILD ANCHOR WATER ARRAY
                    anchorcontacts[anchor]['waters'] = aguas #water atom ids that crossed this anchor
                    anchorcontacts[anchor]['watersite_centroid'] = centroid
            
            print(f"Generating JSON object data for MD simulation {trajectory}")
            res  = {'filename': f'wasimap_{trajectory}.pdb', 'anchor_contacts' : anchorcontacts, 'important_waters': importantwaters}
            
            print("###########################")
            return res




    #Find the residue position where there is a distance jump... for messed soups where everything is inside a single chain
    def find_first_molecule_end_by_peptide_break(self, traj: md.Trajectory, threshold_A: float = 2.0) -> int | None:
        """
        Assumes proteins are in a single chain with increasing residue IDs.
        Returns the residue.index of the last residue of the first protein molecule,
        detected by a break in the peptide bond (C(i) - N(i+1)) distance.
        """
        top = traj.topology
        chain = top.chain(0)

        prot_res = [r for r in chain.residues if r.is_protein]
        if len(prot_res) < 2:
            return prot_res[-1].index if prot_res else None

        xyz = traj.xyz[0]  # first frame (nm)

        for r_i, r_j in zip(prot_res[:-1], prot_res[1:]):
            aC = next((a.index for a in r_i.atoms if a.name == "C"), None)
            aN = next((a.index for a in r_j.atoms if a.name == "N"), None)
            if aC is None or aN is None:
                continue

            dist_nm = np.linalg.norm(xyz[aC] - xyz[aN])
            dist_A = dist_nm * 10.0
            if dist_A > threshold_A:
                return r_i.index

        # No break found: treat as single continuous protein
        return prot_res[-1].index    