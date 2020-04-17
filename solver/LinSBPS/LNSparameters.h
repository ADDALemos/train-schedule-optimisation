#ifndef LNSparameters_h
#define LNSparameters_h

class LNSparameters {
	public:

		LNSparameters() {
			_sbps = false;
			_backjump_shift_strat = false;
			_alternate_backjumps = false;

			_alternate_sbps_ps = false;
			
			_obj_phasing_zero = false;
			_obj_phasing_one = false;
			
			_negative_sbps_chance = 0;

			_sbps_chance = 0;

			_eproc = false;

		}

		bool _eproc;

		bool _sbps;
		bool _backjump_shift_strat;
		bool _alternate_backjumps;

		bool _alternate_sbps_ps;
		bool _obj_phasing_zero;
		bool _obj_phasing_one;

		int _negative_sbps_chance;

		int _sbps_chance;
};

#endif