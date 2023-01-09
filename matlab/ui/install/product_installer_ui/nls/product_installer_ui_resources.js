/* Copyright 2020-2021 The MathWorks, Inc. */
define({
    root: ({
        dialogTitle: 'Matlab Installer',
        installButtonLabel: 'Begin Install',
        downloadButtonLabel: 'Begin Download',
        forwardButtonLabel: 'Next',
        cancelButtonLabel: 'Cancel',
        finishButtonLabel: 'Close',
        guiltForwardButtonLabel: 'I am a student',
        closeLabel: 'OK',
        continueLabel: 'Continue',
        ignoreLabel: 'Ignore',
        licenseAgreementTitle: 'MathWorks License Agreement',
        patentsAndTrademarksTitle: 'Copyrights, Trademarks, and Patents',
        fikWorkflowLabel: 'Install using File Installation Key',
        proxyTitle: 'Proxy Authentication',
        proxyUser: 'Proxy Server User ID:',
        proxyPassword: 'Proxy Server Password:',
        downloadOnlyWorkflowLabel: 'Download Without Installing',
        licenseManagerWorkflowLabel: 'Network License Manager Installation',
        fikLabel: 'Enter File Installation Key',
        fikContextHelp:
            'A File Installation Key identifies the products you can install and allows you to install them without an internet connection.  You must also have your license file for this workflow. You can get both a File Installation Key and a license file from the MathWorks License Center or your license administrator.',
        downloadFolderContextHelp:
            'Download products without installing now so you can install them on one or more computers later.',
        activationKeyContextHelp:
            'An Activation Key links you to another license. You can get an Activation Key from your license administrator.',
        proxyContextHelp:
            'Your machine has been configured to use a proxy server for internet access. Enter the user name and password to continue with a standard online installation. Contact your System Administrator for proxy server credentials.',

        onlineWorkflowButton: 'I want to do a standard install',
        downloadOnlyWorkflowButton: 'I want to download without installing',
        fikWorkflowButton: 'I have a File Installation Key',
        licenseManagerWorkflowButton: 'I want to install network license manager',

        AKToolTipHelpWithAssociation: 'An Activation Key links you to another license.' +
            ' You can get an Activation Key from your license administrator.',
        AKToolTipHelpNoAssociation: 'An Activation Key links you to a license.' +
            ' You can get an Activation Key from your license administrator.',
        downloadOnlyToolTipHelp: 'Download products without installing now so you can install' +
            ' them on one or more computers later.',
        selectLicenseFileForLMToolTipHelp: 'Generate a license file for your server from the' +
            ' MathWorks License Center.',
        selectLicenseFileForToolTipHelp: 'Get the License File from the MathWorks License Center' +
            ' or your license administrator.',
        licenseManagerToolTipHelp: 'This workflow is for License Administrators only. You need' +
            ' to install the Network License manager on a server for Concurrent and Network Named' +
            ' User licenses.',
        akLabel: 'Enter a Activation Key',
        installationFolderLabel: 'Select destination folder',
        downloadFolderLabel: 'Select destination folder',
        platformSelectionLabel: 'Select target platforms',
        browseButtonLabel: 'Browse',
        restoreDefault: 'Restore Default',
        productTableLabel: 'Select products',
        productTableLabelMinimal: 'Select products (recommended products are preselected)',
        productColumnHeader: 'Select All',
        advancedOptionText: 'Advanced Options',
        optionsTitle: 'Select options',
        shortcutLabel: 'Add shortcut to desktop',
        symbolicLinkLabel: 'Create symbolic links to MATLAB scripts in:',
        userExperienceLabel: 'Improve MATLAB by sending user experience information to MathWorks',
        finalLabel: 'Installation Complete',
        finalLabelDownloadOnly: 'Download Complete',
        fikNavTitle: 'File Installation Key',
        folderSelectionNavTitle: 'Installation Folder',
        productSelectionNavTitle: 'Products',
        optionsNavTitle: 'Options',
        confirmationNavTitle: 'Confirmation',
        licenseSelectionTitle: 'Select License File',
        licenseSelectionFolderLabel: '(enter full path to your license file,' +
            ' including file name)',
        entitlementPanelTitle: 'Select license',
        activationPanelTitle: 'Enter Activation Key',
        activationPanelTitleA2AK: 'Enter Activation Key for your license',
        A2AKPanelContextHelp:
            'An Activation Key is required for access to this license. You can get an Activation Key from your license administrator.',
        entitlementTableTitle: 'Licenses:',
        entitlementTableLicenseCol: 'License',
        entitlementTableLabelCol: 'Label',
        entitlementTableOptionCol: 'License Use and Option',
        activationKeyTitle: 'Enter Activation Key:',
        activationModePanelTitle: 'Choose an authentication option',
        activationModeRadioLabel: 'Authorize this computer now',
        activationModeRadioText:
            'Activate the software for this computer. No internet connection required to run the software.',
        activationModeAuthRadioLabel: 'Authenticate (sign in) each time the software starts',
        activationModeAuthRadioText: 'Internet connection required to run the software.',
        configureService: 'Configure as service',
        configureServiceText: 'On Windows systems, the installer configures the license manager' +
            ' to start automatically when you start your system. If you accept this default' +
            ' configuration, the easiest way to start the license manager is to restart the' +
            ' computer on which you installed the license manager.',
        confirmUserTitle: 'Confirm User',
        confirmUserNameLabel: 'Name',
        confirmUserEmailLabel: 'Email',
        confirmWinUserNameLabel: 'Windows User Name',
        editWinUserNameLabel: 'edit',
        ConfirmLinuxUserName: 'Login Name for glnxa64',
        ConfirmMacUserName: 'Login Name for maci',
        linuxOrMacUserNameLabel: 'Login Name',
        firstName: 'First name',
        lastName: 'Last name',
        getUserInfoBtnLabel: 'I will use this software',
        openDownloadFolder: 'Open Downloads folder',
        activateMATLAB: 'Activate your software',
        activateNote: 'Note: To use the software, you must activate it first.',
        provideUserInfoTitle: 'Provide user information',
        provideUserInfoText: 'Only one person can use this license. Specify the Licensed End User:',
        provideUserInfo: '(if necessary a MathWorks Account will be created for this person.)',
        productDependecyDialogTitle: 'Product Dependencies',
        productDependecyDialogAddBtnLabel: 'Add',
        productDependecyDialogDoNotAddBtnLabel: 'Don\'t add',
        productDependecyDialogCancelBtnLabel: 'Cancel',
        student: 'Yes, I am student',
        notStudent: 'No, I am not student',
        studentConfirmTitle: 'Student Use Policy',
        studentConfirmMessage: 'This is a Student License. You must be a student enrolled at a' +
            ' degree-granting higher education institution, or a student or teacher at a high' +
            ' school/pre-university school. It is not for commercial or other use. MathWorks' +
            ' offers this license as a special service to students and asks for help in seeing that' +
            ' its terms are not abused.\n\n' +
            'Are you a student using the software for coursework at a school, college, or' +
            ' university?\n\n',
        yes: 'Yes',
        no: 'No',
        userNameModifyTitle: 'Windows User Name Warning',
        userNameModifyMessage1: 'Change ONLY if the software will be run under a' +
            ' different user name than: ${0}',
        userNameModifyMessage2: ' \n \n' +
            'For more information about user names, see this <a' +
            ' href="https://www.mathworks.com/pi_unc_mpi_${0}_${1}" target=\"_blank\">MATLAB' +
            ' Answer</a>.',
        loginNameModifyTitle: 'Login Name Warning',
        loginNameModifyMessage1: 'Change ONLY if the software will be run under a' +
            ' different login name than: ${0}',
        loginNameModifyMessage2: '\n \n' +
            'For more information about login names, see this <a' +
            ' href="https://www.mathworks.com/pi_unc_mpi_${0}_${1}" target=\"_blank\">MATLAB' +
            ' Answer</a>',
        userNameRequired: 'User name is required',
        confirmationDialogTitle: 'Quit Installer',
        confirmationDialogText: 'If you quit now, you must restart the installation from the' +
            ' beginning.<br><br>Do you really want to quit now?',
        mySelectionButton: 'Use my selection',
        recommendedSelectionButton: 'Use recommendation',
        radioButtonDescriptionForLicenseAgreement:
            'Do you accept the terms of the license agreement?',
        radioButtonsYesText: 'Yes',
        radioButtonsNoText: 'No',
        signInButtonText: 'Sign In',
        connectionErrorTitle: 'Connection Error',
        connectionErrorMsg: 'The application could not connect to MathWorks.' +
            '<br>' +
            'To' +
            ' resolve this issue, contact <a href="http://www.mathworks.com/pi_dlerr_mpi "' +
            ' target="_blank">Customer Support</a>',
        confirmationPanelTitle: 'Confirm selections',
        confirmDestinationFolderLabel: 'DESTINATION',
        confirmPlatformsLabel: 'PLATFORMS',
        confirmProductsLabel: 'PRODUCTS',
        confirmLicensingLabel: 'LICENSING',
        installationNotesTitle:
            '<b>Your installation may require additional configuration steps.</b><br><br>',
        licensingNavTitle: 'Licensing',
        destinationNavTitle: 'Destination',
        platformsNavTitle: 'PLATFORMS',
        genericErrorTitle: 'Something Unexpected Occurred.',
        genericErrorMessage:
            'To resolve this issue, contact  <a href="https://www.mathworks.com/support"' +
            ' target="_blank">Technical Support</a>',
        closeButtonLabel: 'Close',
        mathworksTrademarksText:
            'MATLAB and Simulink are registered trademarks of The MathWorks, Inc. ' +
            'Please see mathworks.com/trademarks for a list of additional trademarks. ' +
            'Other product or brand names may be trademarks or registered trademarks of their respective holders.',
        mathworksPatentsText:
            'MathWorks products are protected by patents (see mathworks.com/patents) and copyright laws. ' +
            'Any unauthorized use, reproduction, or distribution may result in civil and criminal penalties.',
        symbolicLinkFolderLabel: 'Select symbolic link folder',
        noArchiveErrorMsg:
            'Unable to locate the product files for installation. For help, see this <a href="http://www.mathworks.com/pi_miferr_mpi "' +
            ' target="_blank">MATLAB Answer</a>.',
        noArchiveErrorTitle: 'Missing Installation Files',
        noInternetNoArchiveErrorMsg1: 'The issue may be one of the following:',
        noInternetNoArchiveErrorMsg2:
            '<li>The installer is unable to connect to MathWorks. This may be because of a firewall or antivirus software.  For more information, see this <a href="http://www.mathworks.com/pi_icerr_mpi "' +
            ' target="_blank">MATLAB Answer</a>.</li>',
        noInternetNoArchiveErrorMsg3:
            '<li>Product files required for offline installation are missing. For more information, see this <a href="http://www.mathworks.com/pi_miferr_mpi "' +
            ' target="_blank">MATLAB Answer</a>.</li>',
        noInternetNoArchiveErrorTitle: 'Unable to continue'
    }),
    'ja-jp': true,
    'ko-kr': true,
    'zh-cn': true
});
