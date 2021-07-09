def generate_gaussianFile(geom, grid, logger, outdir="./", igrid=0, maxbq=200):
    gaussianfile = outdir + \
        "input_batch_{:05d}.com".format(igrid)
    f = open(gaussianfile, "w")
#    f.write("%OldChk=/home/aartigas/chk/molecule_spe.chk\n".format())
    f.write("%nproc=8\n".format())
    f.write("%mem=1000MB\n".format())
#    f.write("#P b3lyp/6-311++G(d,p) SCF(Tight) CPHF(Separate) Int(Grid=SuperFine) Guess=Read NMR geom=connectivity\n\nTitle\n\n0 1\n".format())
    f.write("#P b3lyp/6-311++G(d,p) SCF(Tight) CPHF(Separate) Int(Grid=SuperFine) NMR geom=connectivity\n\nTitle\n\n0 1\n".format())
    nat = 0
    for at in geom.atoms:
        f.write("{:4s} {:16.10f} {:16.10f} {:16.10f}\n".format(at['label'], at['x'], at['y'], at['z']))
        nat = nat + 1
    nbq = 0
    for at in grid[igrid:]:
        f.write(
            "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f}\n".format(at))
        nbq = nbq + 1
        nat = nat + 1
        igrid = igrid + 1
        if (nbq == maxbq):
            logger.info("Batch generation : {}".format(igrid))
            generate_gaussianFile(
                geom, grid, logger, outdir=outdir, igrid=igrid, maxbq = maxbq)
            break
    f.write("\n")
    for i in range(nat):
        f.write("{}\n".format(i + 1))
    f.write("\n")
    f.close()
    return


