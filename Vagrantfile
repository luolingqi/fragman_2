# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.provision :ansible do |ansible|
    ansible.playbook = "ansible/playbook.yml"
  end

  config.vm.define :web do |web_config|
    # Every Vagrant virtual environment requires a box to build off of.
    web_config.vm.box = "ubuntu/trusty64"

    web_config.vm.host_name = "web"
    web_config.vm.network "private_network", ip: "10.10.10.30"

    web_config.vm.synced_folder "./storage", "/var/storage", create: true
  end

  config.vm.define :db do |db_config|
    db_config.vm.box = "ubuntu/trusty64"

    db_config.vm.host_name = "db"
    db_config.vm.network "private_network", ip: "10.10.10.31"
  end

  config.vm.provider "virtualbox" do |vb|
    # Use VBoxManage to customize the VM. For example to change memory:
    vb.customize ["modifyvm", :id, "--memory", "2048"]
  end
end
