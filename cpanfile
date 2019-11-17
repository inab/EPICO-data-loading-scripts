requires 'Config::IniFiles';
requires 'Net::FTP::AutoReconnect';
requires 'Set::Scalar';
requires 'SQL::Statement';

# These dependencies are from the DarkPAN
requires 'TabParser', '0.01';

requires 'BP::Model', 'v1.1.1';

requires 'BP::Loader', 'v1.0.3';

on test => sub {
    requires 'Test::More', '0.96';
};

on develop => sub {
    requires 'Dist::Milla', '1.0.20';
    requires 'Dist::Zilla::Plugin::MakeMaker';
    requires 'Dist::Zilla::Plugin::Prereqs::FromCPANfile';
    requires 'Dist::Zilla::Plugin::Run', '0.048';
    requires 'OrePAN2';
};
