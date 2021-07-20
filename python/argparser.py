def add_args(parser):
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
    parser.add_argument(
        '--mode',
        type=str,
        default='euler',
        choices=['euler', 'rk4'],
        help='Method used to track the particle'
    )
    return parser.parse_known_args()
