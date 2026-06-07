extern "C" {
#include "lib/mcce.h"
}

#include <stdio.h>
//#include "lib/mcce.h"
//#include <mpi.h>

void welcome();

int main(int argc, char *argv[])
{
   /* Welcome */
   welcome();
/* add */

   /* Do step 0, initialization */
   printf("Step 0. Initialize enviroment\n"); fflush(stdout);
   if (init()) {
      printf("Help message: double check file \"run.prm\" in current directory.\n");
      return USERERR;
   }
   else printf("Step 0 Done.\n\n");


   /* Do step 1, premcce */
   if (env.do_premcce) {
      printf("Step 1. Test and format structral file\n"); fflush(stdout);
      if (premcce()) {return USERERR;}
      else printf("Step 1 Done.\n\n");
   }
   else printf("Not doing \"Step 1. Test and format structral file\"\n\n");


   /* Do step 2. rotamers */
   if (env.do_rotamers) {
      printf("Step 2. Make multi side chain conformers\n"); fflush(stdout);
	if (rotamers()) {
		return USERERR;
    }
    else printf("Step 2 Done.\n\n");
   }
   else printf("Not doing \"Step 2. Make multi side chain conformers\"\n\n");

   /* Do step 3. energies */
   if (env.do_energies) {
      printf("Step 3. Compute energy lookup table\n"); fflush(stdout);
      if (energies()) { return USERERR;}


      else printf("Step 3 Done.\n\n");
   }
   else printf("Not doing \"Step 3. Compute energy lookup table\"\n\n");

   /* Do step 4. Monte Carlo */
   if (env.do_monte) {
      printf("Step 4. Monte Carlo Sampling\n"); fflush(stdout);
      if (!env.monte_adv_opt) {
      if (monte()) {return USERERR;}
           else printf("Step 4 Done.\n\n");
       }
       else {
           if (monte2()) { return USERERR;}
           else printf("Step 4 Done.\n\n");
       }
   }
   else printf("Not doing \"Step 4. Monte Carlo Sampling\"\n\n");


/* Add H bond analysis---Cai */
   if (env.do_analysis){
      printf("Step 6. Hydrogen Bond Network Analysis\n"); fflush(stdout);
      if (analysis_adv()) { return USERERR;}
      else printf("Step 6 Done.\n\n");
   }
   else printf("Not doing \"Step 6. Hydrogen Bond Network Analysis\"\n\n");


   return 0;
}

void welcome()
{  printf("================================================================\n");
   printf("<<< MCCE Multi-Conformation Continuum Electrostatics >>>        \n");
   printf(" Marilyn Gunner's Lab at City College of New York, 2023         \n");
   printf("----------------------------------------------------------------\n");
   printf("Version:         4.0                                            \n");
   printf("MCCE Home:                                                      \n");
   printf("  Private (Dev): https://github.com/GunnerLab/MCCE4             \n");
   printf("  Public:        https://github.com/GunnerLab/MCCE4-Alpha       \n");
   printf("  Docs:          https://gunnerlab.github.io/mcce4_tutorial/    \n");
   printf("Support:         https://github.com/GunnerLab/MCCE4-Alpha/issues\n");
   printf("Developed by:    The Gunner Lab                                 \n");
   printf("                 https://gunnerlab.ccny.cuny.edu/               \n");
   printf("Citations:       If you publish data calculated with MCCE, please cite the papers in\n");
   printf("                 https://gunnerlab.github.io/mcce4_tutorial/docs/mcce/\n");
   printf("================================================================\n");
   fflush(stdout);

   return;
}
