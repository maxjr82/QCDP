
import os
import re
import sys
import argparse
import numpy as np

class Gaussian09:

   def __repr__(self):
        return "Extract information from Gaussian09 outputs"

   def __init__(self,output):
       self.out = open(output, 'r')
       self.lines = self.out.readlines()
       self.out.close()

   @property 
   def get_atom_labels(self):
      out_string = ''.join(self.lines)
      tmp = re.findall(r'^ [A-Z]   ', out_string, flags=re.MULTILINE)
      labels = np.asarray([i.strip() for i in tmp])
      return labels

   @property
   def get_coords(self):
      xyz_list = list()
      match = False
      count = 0
      for line in self.lines:
         if 'Input orientation' in line:
            match = True
            count += 1
            continue
         elif match and (count >= 5):
            if len(line.split()) != 6:
               break
            xyz = np.array(line.split()[-3:], dtype=np.float64)
            xyz = xyz.reshape(-1,3)
            xyz_list.append(xyz) 
         elif match:
            count += 1
      coords = np.vstack(xyz_list)
      coords = np.expand_dims(coords, axis=0)
      return coords

   @property
   def get_forces(self):
      ua_to_kcal = 1185.8212143783846
      start_reading = False
      count = 0
      F_list = list()
      for line in self.lines:
         if start_reading:
            count += 1
         if 'Forces (Hartrees/Bohr)' in line:
            start_reading = True
            count += 1
            continue
         elif start_reading and (count >= 4):
            if len(line.split()) != 5:
               break
            xyz = np.array(line.split()[-3:], dtype=np.float64)
            xyz = xyz.reshape(-1,3)
            F_list.append(xyz)
      forces = np.vstack(F_list)
      forces *= ua_to_kcal
      forces = np.expand_dims(forces, axis=0)
      return forces

   @property
   def get_total_energy(self):
      for line in self.lines:
         if 'SCF Done' in line:
            total_energy = np.array(line.split()[4], dtype=np.float64, ndmin=2)
            total_energy *= 627.5096080305927
      return total_energy

   @property
   def get_homo_lumo(self):
      for line in self.lines:
         if ' occ. ' in line:
            e_homo = float(line.split()[-1])
         if ' virt. ' in line:
            e_lumo = float(line.split()[4])
            break
      e_homo_lumo = np.array([e_homo, e_lumo])
      e_homo_lumo = e_homo_lumo.reshape(-1,2)
      e_homo_lumo *= 27.211399 # convert from Hartree to eV
      return e_homo_lumo

   @property
   def get_dipole(self):
   #    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[D]\ *-?\ *[0-9]+)?')
      scinot = re.compile('[-+]?[\d]+\.?[\d]*[Dd](?:[-+]?[\d]+)?')
      for line in self.lines:
         if 'Dipole    ' in line:
            dipole = [x.replace('D','E') for x in re.findall(scinot, line)]
            dipole = np.array(dipole, dtype=np.float64)
            dipole = dipole.reshape(1,3)
            return dipole

   @property
   def get_mulliken(self):
      count = -1
      start_reading = False
      charges = list()
      pattern1 = 'Mulliken charges:'
      pattern2 = 'Mulliken atomic charges:'
      for line in self.lines:
         if pattern1 in line or pattern2 in line:
            start_reading = True
            continue
         if start_reading:
            count += 1
            if count >= 1:
               if len(line.split()) != 3:
                  break
               charges.append(line.split()[-1])
      charges = np.array(charges, dtype=np.float64)
      n_atoms = len(charges)
      charges = charges.reshape(1, n_atoms, 1)
      return charges

   @property
   def get_zpve(self):
      start_reading = False
      for line in self.lines:
         if 'vibrational energy' in line:
            start_reading = True
            continue
         if start_reading:
            zpve = np.array(line.split()[0], dtype=np.float64, ndmin=2)
            return zpve

   @property
   def _get_exact_polar(self):
      for line in self.lines:
         if 'Exact polarizability' in line:
            polar = np.array(line.split()[-6:], dtype=np.float64)
            polar = polar.reshape(-1,6)
            return polar

   @property
   def get_polarizability(self):
      start_reading = False
      scinot = re.compile('[-+]?[\d]+\.?[\d]*[Dd](?:[-+]?[\d]+)?')
      for line in self.lines:
         if 'Polarizability' in line:
            start_reading = True
            p1 = [x.replace('D','E') for x in re.findall(scinot, line)]
         if start_reading:
            p2 = [x.replace('D','E') for x in re.findall(scinot, line)]
            polar = np.array((p1 + p2), dtype=np.float64)
            polar = polar.reshape(-1,6)
            return polar

   @property
   def get_rot_constants(self):
      for line in self.lines:
         if 'Rotational constants' in line:
            rot_const = np.array(line.split()[-3:], dtype=np.float64)
            rot_const = rot_const.reshape(1,3)
            return rot_const

   @property
   def get_elec_spatial_ext(self):
      for line in self.lines:
         if 'Electronic spatial extent' in line:
            r2 = np.array(line.split()[-1], dtype=np.float64, ndmin=2)
            return r2

   @property
   def get_thermal_energies(self):
      count = 0 
      thermal_energies = list()
      pattern = 'Sum of electronic and thermal'
      for line in self.lines:
         if pattern in line:
            thermal_energies.append(line.split()[-1])
            count += 1
         if count == 3:
            break
      thermal_energies = np.array(thermal_energies, dtype=np.float64)
      thermal_energies *= 627.5096080305927
      thermal_energies = thermal_energies.reshape(-1,3)
      return thermal_energies

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-log", default=None, type=str, help = "full path to the output file")
    parser.add_argument("-prop", default=None, type=str, help = "property to be read from the output")

    if len(sys.argv) == 1:
        args = parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    log_file = args.log
    property = args.prop

    try:
       g09_output = Gaussian09(log_file)
       print(g09_output)
       print(" ")
    except:
       print("Oops! File not found. Try again.")
       sys.exit()

    all_properties = {'labels': g09_output.get_atom_labels,
                      'coords': g09_output.get_coords,
                      'Etot': g09_output.get_total_energy,
                      'e_hl': g09_output.get_homo_lumo,
                      'alpha': g09_output.get_polarizability,
                      'dipole': g09_output.get_dipole,
                      'e_vib': g09_output.get_zpve,
                      'rot': g09_output.get_rot_constants,
                      'r2': g09_output.get_elec_spatial_ext,
                      'forces': g09_output.get_forces,
                      'e_thermal': g09_output.get_thermal_energies,
                      'mulliken': g09_output.get_mulliken}

    if property is None or property == 'all':
        print('Reading all properties from file: {}'.format(log_file))
        print('-------------------------------------')
        for k, val in all_properties.items():
           print(k, val)
    else:
       if ',' in property:
          props = [p.strip() for p in property.split(',')]
          print('Reading {} properties from file: {}'.format(len(props),log_file))
          print('-------------------------------------')
          for p in props:
             try:
                print(p, all_properties[p])
             except:
                print(" ")
                print("Oops! Key not found in dictionary.\n")
                print("The available options are", list(all_properties.keys()))
       else:
           p = property.strip()
           print('Reading {} from file: {}'.format(p,log_file))
           print('-------------------------------')
           try:
               print(p, all_properties[p])
           except:
               print(" ")
               print("Oops! Key not found in dictionary.\n")
               print("The available options are", list(all_properties.keys()))
