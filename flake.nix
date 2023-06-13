{
  description = "simple python environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils }: (utils.lib.eachSystem ["x86_64-linux" ] (system:
    let
      pkgs = nixpkgs.legacyPackages.${system};
      pypkgs = pkgs.python3Packages;
      pyenv = pkgs.mkShell rec {
        pname = "simple";
        nativeBuildInputs = with pypkgs; [
          scipy
          numpy
          matplotlib
          tkinter
          pandas
          jupyterlab
          notebook
        ];
        shellHook = ''
          # export NIX_SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
          # export OMPI_MCA_plm_rsh_agent=${pkgs.openssh}/bin/ssh

          SOURCE_DATE_EPOCH=$(date +%s)
          export PYTHONUSERBASE=$PWD/.local
          export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
          export PYTHONPATH=$PYTHONPATH:$USER_SITE
          export PATH=$PATH:$PYTHONUSERBASE/bin

          export TK_LIBRARY="${pkgs.tk}/lib/${pkgs.tk.libPrefix}"
        '';
      };
    in
      {
        devShells.default = pyenv;
      }
  ));
}
