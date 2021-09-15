def add_args(parser, mode):
    parser.add_argument(
        '--x',
        type=float,
        default=0.00,
        help='Initial beam X position [cm].'
    )
    parser.add_argument(
        '--y',
        type=float,
        default=0.00,
        help='Initial beam Y position [cm].'
    )
    parser.add_argument(
        '--energy',
        type=float,
        default=1380.00,
        help='Beam energy [GeV].'
    )
    
    if mode == 'sim':
        parser.add_argument(
            '--step_size',
            type=float,
            default=1380.00,
            help='Step size [cm].'
        )
        parser.add_argument(
            '--mass_interaction',
            type=float,
            default=0.139,
            help='Mass of the interaction TLorentzVector [GeV].'
        )
        parser.add_argument(
            '--npartons',
            type=int,
            default=1,
            help='Number of partons (to scale the energy)',
        )

        parser.add_argument(
            '--mode',
            type=str,
            default='euler',
            choices=['euler', 'rk4'],
            help='Method used to track the particle'
        )

    if mode == 'gen':
        parser.add_argument(
            '--dataset',
            type=str,
            choices=['boltz', 'gaus'],
            help='Distribution to inspect.'
        )
        
    return parser.parse_known_args()
