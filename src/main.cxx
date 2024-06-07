#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <THashList.h>
#include <TList.h>

using namespace std;

#define ENDC "\033[0m"
#define CBOLDRED "\033[1m\033[31m"
#define CBLODYELLOW "\033[1m\033[33m"
#define CRED "\033[31m"
#define CYELLOW "\033[33m"

class MyParticlePDG : public TParticlePDG {
	public:
		MyParticlePDG(const char *name, const char *title, Double_t m, Bool_t stable, Double_t width, Double_t charge, const char *particleClass, Int_t pdgCode, Int_t anti, Int_t trackingCode)
						: TParticlePDG(name, title, m, stable, width, charge, particleClass, pdgCode, anti, trackingCode) {}
		MyParticlePDG(const TParticlePDG &particle)
						: TParticlePDG(particle) {}

		static MyParticlePDG *Create(const TParticlePDG *particle)
		{
			if (!particle)
				return nullptr;
			return new MyParticlePDG(*particle);
		}

		virtual void Print(string opts = "bm") const
		{
			int flg_basic = 0; // name and pid
			int flg_mass = 0;
			int flg_decay = 0;
			if(opts.find("D")!=string::npos || opts.find("d")!=string::npos)
				flg_decay = 1;
			if(opts.find("b")!=string::npos || opts.find("B")!=string::npos)
				flg_basic = 1;
			if(opts.find("m")!=string::npos || opts.find("M")!=string::npos)
				flg_mass = 1;
			if(flg_basic != 0)
				printf("%-20s  %6d\t\n",GetName(),fPdgCode);
			if(flg_mass != 0)
			{
				if (!fStable) {
					printf("Mass:%9.4f Width (GeV):%11.4e\tCharge: %5.1f\n",
							fMass, fWidth, fCharge);
				} else {
					printf("Mass:%9.4f Width (GeV): Stable\tCharge: %5.1f\n",
							fMass, fCharge);
				}
			}
			if (fDecayList && flg_decay != 0) {
				int banner_printed = 0;
				TIter next(fDecayList);
				TDecayChannel* dc;
				while ((dc = (TDecayChannel*)next())) {
					if (! banner_printed) {
						PrintDecayChannel(dc,"banner");
						banner_printed = 1;
					}
					PrintDecayChannel(dc,"data");
				}
			}
		}
};

int LevenshteinDistance(const std::string &s1, const std::string &s2) {
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

    d[0][0] = 0;
    for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for (unsigned int i = 1; i <= len1; ++i)
        for (unsigned int j = 1; j <= len2; ++j)
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1)});
    return d[len1][len2];
}

std::vector<TParticlePDG*> FuzzySearch(TDatabasePDG *db, const std::string &query) {
    std::vector<TParticlePDG*> results;
				db->GetParticle(1);
				const THashList *plist=db->ParticleList();
    TIter next(plist);
    TParticlePDG *particle;
				while ((particle = (TParticlePDG *)next()))
				{
					std::string particleName = particle->GetName();
					int maxDistance = particleName.length() - query.length();
					if (maxDistance < 0)
						maxDistance = 0;
					maxDistance += 2; // allow for up to 3 mismatches
					if(maxDistance > 4)
						maxDistance = 4; // limit the max distance to 5
					int distance = LevenshteinDistance(query, particleName);
					if (distance <= maxDistance)
					{
						results.push_back(particle);
        }
				}
				return results;
}

void helps()
{
 cout << "Usage: \n"
						<< endl;
 cout << "PDGref [-n <particle name>] [-id <particle PDG ID>] [-D] [-m] [-b]\n"
						<< endl;
 cout << "   -n or --name: print the particle name and PDG ID of the particle whose name is similar to the input name.\n"
						<< endl;
 cout << "   -id or --PDGid: print the particle name and PDG ID of the particle whose PDG ID is equal to the input PDG ID.\n"
						<< endl;
 cout << "   -D or --Decay: print the decay channels of the particle, or only print basic information of the particle if -D is not specified.\n"
						<< endl;
	cout << "   -m or --Mass: print the mass and width of the particle, or only print basic information of the particle if -m is not specified.\n"
						<< endl;
 cout << "   -b or --basic: print only name and PDG ID of the particle, without any decay channels or mass information.\n"
						<< "   If there is no -b, -m or -D option, -b -m is the default option.\n "
						<< endl;
 cout << "   -h or --help: print the usage of the program.\n"
						<< endl;
}

int main(int argc, char *argv[])
{
 int flg_help = 0;

 // for input information
	int particle_name_loc  = -1;
	int particle_pdgid_loc = -1;

	// for print information
	int flg_basic_print    = 0;
	int flg_Decay_print    = 0;
	int flg_mass_print    	= 0;

 cout << endl;
 cout << "PDGref v1.0\n"
						<< endl;

 for (int i = 1; i < argc; i++)
 {
		string args = argv[i];
		if (args=="-h"||args=="--help" || argc < 2)
			flg_help=1;
		else if((args == "-n" || args == "--name") && particle_name_loc == -1)
		{
			i++;
			if (i >= argc)
			{
				cerr << CBOLDRED << "Error: " << ENDC << "No particle name spicified for -n/--name option\n"
									<< endl;
				return 1;
			}
			particle_name_loc = i;
		}
		else if ((args == "-id" || args == "--PDGid") && particle_pdgid_loc == -1)
		{
			i++;
			if (i >= argc)
			{
				cerr << CBOLDRED << "Error: " << ENDC << "No PDGid spicified for -id/--PDGid option\n"
									<< endl;
				return 1;
			}
			particle_pdgid_loc = i;
		}
		else if(args == "-D" || args == "--Decay")
		{
			flg_Decay_print=1;
		}
		else if(args == "-m" || args == "--Mass")
		{
			flg_mass_print = 1;
		}
		else if(args == "-b" || args == "--basic")
		{
			flg_basic_print = 1;
		}
		else
		{
			cout << "Invalid argument: " << args << "\n"
								<< endl;
			flg_help=1;
			break;
		}
	}

 if(flg_help == 1)
 {
		helps();
		return 1;
 }

 if(flg_Decay_print == 0 && flg_mass_print == 0 && flg_basic_print == 0)
 {
		flg_basic_print = 1;
		flg_mass_print = 1;
 }

 string opts = "";
 if (flg_basic_print == 1)
		opts += "b";
 if (flg_mass_print == 1)
		opts += "m";
 if (flg_Decay_print == 1)
		opts += "d";

 TDatabasePDG *db = TDatabasePDG::Instance();
 MyParticlePDG *particle;
 string particlename;
 int pdgid;
 if (particle_pdgid_loc != -1)
		pdgid = atoi(argv[particle_pdgid_loc]);
 if (particle_name_loc != -1)
		particlename = argv[particle_name_loc];

 if (particle_pdgid_loc != -1)
 {
		//TParticlePDG *particle_r = db->GetParticle(pdgid);
			particle=MyParticlePDG::Create(db->GetParticle(pdgid));
			if (particle == nullptr)
			{
				cout << CBOLDRED << "Error: " << ENDC << "Invalid PDG ID: " << pdgid << "\n"
									<< endl;
				if (particle_name_loc == -1)
					return 1;
				else
					cout << CBLODYELLOW << "Warning: " << ENDC << "using name search\n"
										<< endl;
		}
		else
		{
			if (particle_name_loc != -1 && particle->GetName() != particlename)
				cout << CBLODYELLOW << "Waring: " << ENDC << "PDG ID doesn't match name, only print PDG ID searched result\n"
									<< endl;

			particle->Print(opts);
			cout << endl;
			return 0;
		}
 }
	// following code need to debug 

 if(particle_name_loc != -1)
	{
		particle = MyParticlePDG::Create(db->GetParticle(particlename.c_str()));
		//particle = new MyParticlePDG(*db->GetParticle(particlename.c_str()));
		if(particle != NULL)
		{
			particle->Print(opts);
			cout << endl;
			return 0;
		}
		else
		{
			cout << CBLODYELLOW << "Warning: " << ENDC << "there is no particle named " << particlename << "\n"
								<< "Searcing for similar name ... \n"
								<< endl;

			vector<TParticlePDG *> sim_particles = FuzzySearch(db, particlename);
			if (sim_particles.size() == 0)
			{
				cerr << CBOLDRED << "Error: " << ENDC << "no similar particle found\n"
									<< endl;
				return 1;
			}
			for (int ip = 0; ip < sim_particles.size(); ip++)
			{
				particle = new MyParticlePDG(*sim_particles[ip]);
				cout << "///////////////////" << endl;
				particle->Print();
				cout << endl;
			}
			return 0;
		}
	}

 cerr << CBOLDRED << "Error: " << ENDC << "There must be something wrong!\n"
						<< endl;
 helps();

 return 0;
}
